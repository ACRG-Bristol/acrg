import functools
import json
from pathlib import Path
from typing import Any, Optional

import xarray as xr

from attribute_parsers import add_variable_attrs, get_data_var_attrs, make_global_attrs
from countries import Countries
from file_processing import get_netcdf_files
from helpers import (
    calc_mode,
    convert_time_to_unix_epoch,
    make_quantiles,
    sparse_xr_dot,
)
from process_rhime_output import InversionOutput


PARIS_FORMATTING_PATH = Path(__file__).parent


def calculate_stats(
    ds: xr.Dataset, name: str, chunk_dim: str, chunk_size: int = 10, var_names: Optional[list[str]]=None
) -> list[xr.Dataset]:
    output = []
    if var_names is None:
        var_names = list(ds.data_vars)
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
    with open(PARIS_FORMATTING_PATH / "iso3166.json", "r") as f:
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


def main(species: str, output_file_path: str, country_files_root: str, min_model_error: float = 0.1):
    files = get_netcdf_files(output_file_path, filename_search=species.upper())
    inv_outs = [InversionOutput.from_rhime(xr.open_dataset(file), min_model_error) for file in files]

    for inv_out in inv_outs:
        inv_out.sample_predictive_distributions()

    # calculate country stats
    country_files_path = Path(country_files_root)
    countries_path = country_files_path / "country_EUROPE.nc"
    countries_ukmo_path = country_files_path / "country-ukmo_EUROPE.nc"

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

    countries = Countries(xr.open_dataset(countries_path), paris_countries)
    countries.merge(Countries(xr.open_dataset(countries_ukmo_path), paris_countries_ukmo))

    country_traces = [
        countries.get_country_trace("sf6", inv_out).expand_dims({"time": [inv_out.start_time()]})
        for inv_out in inv_outs
    ]

    country_traces_merged = xr.concat(country_traces, dim="time")

    country_output = xr.merge(calculate_stats(country_traces_merged, "country", chunk_dim="country", chunk_size=1))

    # calculate flux stats
    traces = [inv_out.get_trace_dataset(convert_nmeasure=False, var_names="x") for inv_out in inv_outs]
    stats = [xr.merge(calculate_stats(trace, "flux", chunk_dim="nx")) for trace in traces]

    flux_stats = [
        sparse_xr_dot((inv_out.flux * inv_out.basis), stats_ds).expand_dims({"time": [inv_out.start_time()]})
        for inv_out, stats_ds in zip(inv_outs, stats)
    ]

    flux_all_times = xr.concat(flux_stats, dim="time")

    flux_attrs = get_data_var_attrs(str(PARIS_FORMATTING_PATH / "netcdf_template_emissions_bm_edits.txt"))

    # apply `get_country_code` to each element of `country` coordinate
    country_codes = list(map(get_country_code, map(str, country_output.country.values)))
    country_output = country_output.assign_coords(country=country_codes)

    emissions = (
        xr.merge([flux_all_times, country_output])
        .pipe(convert_time_to_unix_epoch)
        .rename(probs="quantile")
        .pipe(add_variable_attrs, flux_attrs)
    )

    emissions.attrs = make_global_attrs("flux")

    return emissions
