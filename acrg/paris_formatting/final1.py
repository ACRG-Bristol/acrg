import functools
import json
from pathlib import Path
from typing import Any, cast, Optional, Union

import xarray as xr
from attribute_parsers import add_variable_attrs, get_data_var_attrs, make_global_attrs
from file_processing import get_netcdf_files, get_rhime_outs
from helpers import (
    calc_mode,
    convert_time_to_unix_epoch,
    get_area_grid,
    get_country_trace,
    get_x_to_country_mat,
    get_xr_dummies,
    make_quantiles,
    sparse_xr_dot,
)
from sampling import (
    convert_idata_to_dataset,
    get_prior_samples,
    get_posterior_predictive_samples,
    get_rhime_model,
    get_sampling_kwargs_from_rhime_outs,
    make_idata_from_rhime_outs,
)


def get_all_traces(
    outs: list[xr.Dataset], min_model_error: float, thin_by: Optional[int] = None
) -> list[xr.Dataset]:
    # reconstruct RHIME model for each RHIME output dataset
    model_kwargs = [get_sampling_kwargs_from_rhime_outs(ds, min_model_error=min_model_error) for ds in outs]
    models = [get_rhime_model(ds, **kwargs) for ds, kwargs in zip(outs, model_kwargs)]  # type: ignore

    # reconstruct posterior InferenceData from RHIME output datasets
    idatas = [make_idata_from_rhime_outs(ds) for ds in outs]

    # get prior and prior predictive samples
    for idata, model in zip(idatas, models):
        idata.extend(get_prior_samples(model))

    # get posterior predictive samples
    for idata, model in zip(idatas, models):
        idata.extend(get_posterior_predictive_samples(idata=idata, model=model, thin_by=thin_by))

    result = [convert_idata_to_dataset(idata) for idata in idatas]

    return result


def make_country_traces(
    species: str,
    outs: list[xr.Dataset],
    all_traces: list[xr.Dataset],
    country_files: list[Union[str, Path]],
    country_selections: Optional[list[Optional[list[str]]]] = None,
) -> xr.Dataset:
    # make sure country_selections is the same length as country_files
    if country_selections is not None and len(country_files) != len(country_selections):
        raise ValueError(
            f"List of country selections must have the same length as the list of country files: {len(country_files)} != {len(country_selections)}."
        )
    if country_selections is None:
        _country_selections = [None] * len(country_files)
    else:
        _country_selections = cast(list[Optional[list[str]]], country_selections)

    area_grid = get_area_grid(outs[0].lat, outs[0].lon)

    country_datasets = [xr.open_dataset(file) for file in country_files]

    # filter country datasets by country_selections
    all_country_selections = []
    for country_selection in _country_selections:
        if country_selection is not None:
            all_country_selections.extend(country_selection)

    if len(set(all_country_selections)) < len(all_country_selections):
        raise ValueError(
            "Duplicate countries selected. Make sure lists in `country_selections` are disjoint."
        )

    # get matrices mapping x traces to countries
    all_x_to_country_mats = []
    for country_ds, country_selection in zip(country_datasets, _country_selections):
        x_to_country_mats = [
            get_x_to_country_mat(
                country_ds,
                hbmcmc_outs=outs_ds,
                area_grid=area_grid,
                basis_cat_dim="nparam",
                country_selection=country_selection,
            )
            for outs_ds in outs
        ]
        all_x_to_country_mats.append(x_to_country_mats)

    x_to_country_mats_stacked = [
        xr.concat(mats_tup, dim="ncountries") for mats_tup in zip(*all_x_to_country_mats)
    ]

    # get country traces by multiplying x traces by x-to-country matrices
    country_traces = [
        get_country_trace(species, x_trace=traces[["xtrace_posterior", "x_prior"]], x_to_country=mat)
        for traces, mat in zip(all_traces, x_to_country_mats_stacked)
    ]

    # add times and concatenate
    times = [ds.Ytime.min().values for ds in outs]

    all_country_traces = xr.concat(
        [trace.expand_dims({"time": [time]}) for trace, time in zip(country_traces, times)], dim="time"
    )

    return cast(xr.Dataset, all_country_traces)


def calculate_stats(
    ds: xr.Dataset, name: str, chunk_dim: str, chunk_size: int = 10, var_names=["xtrace_posterior", "x_prior"]
) -> list[xr.Dataset]:
    output = []
    for var_name in var_names:
        suffix = "apost" if "posterior" in var_name else "apriori"
        stats = [
            ds[var_name].mean("draw").rename(f"{name}{suffix}"),
            calc_mode(ds[var_name].dropna(dim="draw").chunk({chunk_dim: chunk_size}), sample_dim="draw")
            .compute()
            .rename(f"{name}{suffix}_mode"),
            make_quantiles(ds[var_name].dropna(dim="draw"), sample_dim="draw").rename(f"q{name}{suffix}"),
        ]
        output.extend(stats)
    return output


@functools.lru_cache
def get_iso3166_codes() -> dict[str, Any]:
    with open("/home/brendan/Documents/acrg/acrg/paris_formatting/iso3166.json", "r") as f:
        iso3166 = json.load(f)
    return iso3166


def get_country_code(x: str, iso3166: Optional[dict[str, dict[str, Any]]] = None) -> str:
    if iso3166 is None:
        iso3166 = get_iso3166_codes()

    for k, v in iso3166.items():  # type: ignore
        names = [v["iso_long_name"].lower()] + [name.lower() for name in v["unofficial_names"]]
        if any(x.lower() in name for name in names):
            return k

    return x


def main():
    species = "sf6"
    files = get_netcdf_files("/home/brendan/Documents/inversions/plotting/sf6_best", filename_search="SF6")
    outs = get_rhime_outs(files)[:2]  # limit to 2 files for testing

    all_traces = get_all_traces(outs, min_model_error=0.1, thin_by=10)

    inversions_path = Path("/home/brendan/Documents/inversions/openghg_inversions/")
    countries_path = inversions_path / "countries" / "country_EUROPE.nc"
    countries_ukmo_path = inversions_path / "countries" / "country-ukmo_EUROPE.nc"

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

    country_traces_merged = make_country_traces(
        species,
        outs,
        all_traces,
        [countries_path, countries_ukmo_path],
        [paris_countries, paris_countries_ukmo],
    )

    country_output = xr.merge(calculate_stats(country_traces_merged, "country", "ncountries", 1))

    fluxes = [ds.fluxapriori for ds in outs]
    basis_mats = [get_xr_dummies(ds.basisfunctions, cat_dim="nparam") for ds in outs]
    traces = [trace[["xtrace_posterior", "x_prior"]] for trace in all_traces]

    stats = [xr.merge(calculate_stats(trace, "flux", "nparam")) for trace in traces]

    flux_stats = [
        sparse_xr_dot((flux * mat), stats_ds) for flux, mat, stats_ds in zip(fluxes, basis_mats, stats)
    ]

    times = [ds.Ytime.min().values for ds in outs]
    flux_stats_plus_time = [fs.expand_dims({"time": [time]}) for fs, time in zip(flux_stats, times)]
    flux_all_times = xr.concat(flux_stats_plus_time, dim="time")

    flux_attrs = get_data_var_attrs(
        "/home/brendan/Documents/acrg/acrg/paris_formatting/netcdf_template_emissions_bm_edits.txt"
    )

    # apply `get_country_code` to each element of `country` coordinate
    country_codes = list(map(get_country_code, map(str, country_output.country.values)))
    country_output = country_output.swap_dims(ncountries="country").assign_coords(country=country_codes)

    emissions1 = xr.merge([flux_all_times, country_output])

    emissions2 = convert_time_to_unix_epoch(emissions1)
    emissions3 = emissions2.rename(probs="quantile")

    emissions4 = add_variable_attrs(emissions3, flux_attrs)

    emissions4.attrs = make_global_attrs("flux")

    return emissions4
