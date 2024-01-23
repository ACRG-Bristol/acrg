import json
from pathlib import Path
from typing import Any, Optional, Union

import xarray as xr
from attribute_parsers import add_variable_attrs, get_data_var_attrs, make_global_attrs
from file_processing import get_netcdf_files, get_rhime_outs
from helpers import (
    calc_mode,
    convert_time_to_unix_epoch,
    get_area_grid_data_array,
    get_country_trace,
    get_prior_samples,
    get_x_to_country_mat,
    get_xr_dummies,
    make_quantiles,
    sparse_xr_dot,
)


species = "sf6"
files = get_netcdf_files("/home/brendan/Documents/inversions/plotting/sf6_best", filename_search="SF6")

outs = get_rhime_outs(files)


def get_prior_from_attrs(x: str) -> dict[str, Union[str, float]]:
    """Get dict of prior arguments from RHIME output attributes."""
    prior = {}
    split = x.split(",")
    for k, v in zip(split[::2], split[1::2]):
        prior[k] = v if k == "pdf" else float(v)
    return prior


sampling_kwargs = dict(
    xprior=get_prior_from_attrs(outs[0].attrs["Emissions Prior"]),
    bcprior=get_prior_from_attrs(outs[0].attrs["BCs Prior"]),
    min_model_error=0.15,
)

prior_samples = [get_prior_samples(ds, **sampling_kwargs).prior.squeeze(drop=True) for ds in outs]


def combine_outs_and_prior_samples(outs: xr.Dataset, prior_samples: xr.Dataset) -> xr.Dataset:
    """Combine trace variables from RHIME outputs and prior samples."""
    trace_dvs = [dv for dv in outs.data_vars if "trace" in str(dv)]
    prior_rename_dict = {dv: str(dv) + "_prior" for dv in prior_samples.data_vars}
    return xr.merge([outs[trace_dvs], prior_samples.rename_vars(prior_rename_dict)])


all_traces = [
    combine_outs_and_prior_samples(outs_ds, prior_samples_ds)
    for outs_ds, prior_samples_ds in zip(outs, prior_samples)
]

inversions_path = Path("/home/brendan/Documents/inversions/openghg_inversions/")
countries = xr.open_dataset(inversions_path / "countries" / "country_EUROPE.nc")
countries_ukmo = xr.open_dataset(inversions_path / "countries" / "country-ukmo_EUROPE.nc")


area_grid = get_area_grid_data_array(outs[0].lat, outs[0].lon)

country_mats = [
    get_x_to_country_mat(countries, hbmcmc_outs=outs_ds, area_grid=area_grid, basis_cat_dim="nparam")
    for outs_ds in outs
]
country_ukmo_mats = [
    get_x_to_country_mat(countries_ukmo, hbmcmc_outs=outs_ds, area_grid=area_grid, basis_cat_dim="nparam")
    for outs_ds in outs
]


country_traces = [
    get_country_trace("sf6", x_trace=traces[["xtrace", "x_prior"]], x_to_country=mat)
    for traces, mat in zip(all_traces, country_mats)
]
country_ukmo_traces = [
    get_country_trace("sf6", x_trace=traces[["xtrace", "x_prior"]], x_to_country=mat)
    for traces, mat in zip(all_traces, country_ukmo_mats)
]

paris_countries = [
    "BELGIUM",
    "SWITZERLAND",
    "AUSTRIA",
    "ITALY",
    "NETHERLANDS",
    "CZECHIA",
    "POLAND",
    "HUNGARY",
    "SLOVAKIA",
    "SWEDEN",
    "FINLAND",
]
country_filt = countries.name.isin(paris_countries)

paris_countries_ukmo = [
    "BENELUX",
    "RestEU",
    "SpaPor",
    "IRELAND",
    "UNITED KINGDOM",
    "FRANCE",
    "GERMANY",
    "DENMARK",
    "NORWAY",
]
country_ukmo_filt = countries_ukmo.name.isin(paris_countries_ukmo)

times = [ds.Ytime.min().values for ds in outs]

country_traces_concat = xr.concat(
    [
        trace.where(country_filt, drop=True).expand_dims({"time": [time]})
        for trace, time in zip(country_traces, times)
    ],
    dim="time",
)
country_ukmo_traces_concat = xr.concat(
    [
        trace.where(country_ukmo_filt, drop=True).expand_dims({"time": [time]})
        for trace, time in zip(country_ukmo_traces, times)
    ],
    dim="time",
)
country_traces_merged = xr.concat([country_traces_concat, country_ukmo_traces_concat], dim="ncountries")


def calculate_stats(ds: xr.Dataset, name: str, chunk_dim: str, chunk_size: int = 10) -> list[xr.Dataset]:
    output = [
        ds.xtrace.mean("draw").rename(f"{name}apost"),
        calc_mode(ds.xtrace.chunk({chunk_dim: chunk_size}), sample_dim="draw")
        .compute()
        .rename(f"{name}apost_mode"),
        make_quantiles(ds.xtrace, sample_dim="draw").rename(f"q{name}apost"),
        ds.x_prior.mean("draw").rename(f"{name}apriori"),
        calc_mode(ds.x_prior.dropna(dim="draw").chunk({chunk_dim: chunk_size}), sample_dim="draw")
        .compute()
        .rename(f"{name}apriori_mode"),
        make_quantiles(ds.x_prior.dropna(dim="draw"), sample_dim="draw").rename(f"q{name}apriori"),
    ]
    return output


country_output = xr.merge(calculate_stats(country_traces_merged, "country", "ncountries", 1))


fluxes = [ds.fluxapriori for ds in outs]
basis_mats = [get_xr_dummies(ds.basisfunctions, cat_dim="nparam") for ds in outs]
traces = [trace[["xtrace", "x_prior"]] for trace in all_traces]

stats = [xr.merge(calculate_stats(trace, "flux", "nparam")) for trace in traces]

flux_stats = [sparse_xr_dot((flux * mat), stats_ds) for flux, mat, stats_ds in zip(fluxes, basis_mats, stats)]

flux_stats_plus_time = [fs.expand_dims({"time": [time]}) for fs, time in zip(flux_stats, times)]
flux_all_times = xr.concat(flux_stats_plus_time, dim="time")


flux_attrs = get_data_var_attrs(
    "/home/brendan/Documents/acrg/acrg/paris_formatting/netcdf_template_emissions_bm_edits.txt"
)


def get_country_code(x: str, iso3166: Optional[dict[str, dict[str, Any]]] = None) -> str:
    if iso3166 is None:
        with open("iso3166.json", "r") as f:
            iso3166 = json.load(f)

    for k, v in iso3166.items():  # type: ignore
        names = [v["iso_long_name"].lower()] + [name.lower() for name in v["unofficial_names"]]
        if any(x.lower() in name for name in names):
            return k

    return x


with open("/home/brendan/Documents/acrg/acrg/paris_formatting/iso3166.json", "r") as f:
    iso3166 = json.load(f)

country_output = country_output.swap_dims(ncountries="country").assign_coords(
    country=list(map(lambda x: get_country_code(x, iso3166), map(str, country_output.country.values)))
)

emissions1 = xr.merge([flux_all_times, country_output])

emissions2 = convert_time_to_unix_epoch(emissions1)
emissions3 = emissions2.rename_vars(probs="quantile").swap_dims(probs="quantile")


emissions4 = add_variable_attrs(emissions3, flux_attrs)

emissions4.attrs = make_global_attrs("flux")
