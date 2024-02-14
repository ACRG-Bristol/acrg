import functools
import json
from pathlib import Path
from typing import Any, Optional

import xarray as xr

from attribute_parsers import add_variable_attrs, get_data_var_attrs, make_global_attrs
from countries import Countries
from file_processing import get_netcdf_files
from helpers import calc_mode, convert_time_to_unix_epoch, make_quantiles, sparse_xr_dot
from process_rhime_output import InversionOutput

PARIS_FORMATTING_PATH = Path(__file__).parent


def calculate_stats(
    ds: xr.Dataset,
    name: str,
    chunk_dim: str,
    chunk_size: int = 10,
    var_names: Optional[list[str]] = None,
    report_mode: bool = False,
    add_bc_suffix: bool = False,
) -> list[xr.Dataset]:
    output = []
    if var_names is None:
        var_names = list(ds.data_vars)
    for var_name in var_names:
        suffix = "apost" if "posterior" in var_name else "apriori"
        if add_bc_suffix:
            suffix += "BC"

        if report_mode:
            stats = [
                calc_mode(ds[var_name].dropna(dim="draw").chunk({chunk_dim: chunk_size}), sample_dim="draw")
                .compute()
                .rename(f"{name}{suffix}"),
            ]
        else:
            stats = [
                ds[var_name].mean("draw").rename(f"{name}{suffix}"),
                calc_mode(ds[var_name].dropna(dim="draw").chunk({chunk_dim: chunk_size}), sample_dim="draw")
                .compute()
                .rename(f"{name}{suffix}_mode"),
            ]
        stats.append(
            make_quantiles(ds[var_name].dropna(dim="draw"), sample_dim="draw").rename(f"q{name}{suffix}")
        )

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


def get_inversion_outputs_with_samples(
    species: str, output_file_path: str, min_model_error: float = 0.1
) -> list[InversionOutput]:
    """Create a list of InversionOutputs given a path to RHIME inversion outputs."""
    files = get_netcdf_files(output_file_path, filename_search=species.upper())
    inv_outs = [InversionOutput.from_rhime(xr.open_dataset(file), min_model_error) for file in files]

    for inv_out in inv_outs:
        inv_out.sample_predictive_distributions()

    return inv_outs


def make_country_output(inv_outs: list[InversionOutput], country_files_root: str) -> xr.Dataset:
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

    country_output = xr.merge(
        calculate_stats(country_traces_merged, "country", chunk_dim="country", chunk_size=1)
    )

    # apply `get_country_code` to each element of `country` coordinate
    country_codes = list(map(get_country_code, map(str, country_output.country.values)))
    country_output = country_output.assign_coords(country=country_codes)

    return country_output


def make_flux_outputs(inv_outs: list[InversionOutput]) -> xr.Dataset:
    """Make flux output dataset"""

    # calculate stats on flux traces
    traces = [inv_out.get_trace_dataset(convert_nmeasure=False, var_names="x") for inv_out in inv_outs]
    stats = [xr.merge(calculate_stats(trace, "flux", chunk_dim="nx")) for trace in traces]

    # multiply stats by matrix mapping basis regions to lat/lon
    flux_stats = [
        sparse_xr_dot((inv_out.flux * inv_out.basis), stats_ds).expand_dims({"time": [inv_out.start_time()]})
        for inv_out, stats_ds in zip(inv_outs, stats)
    ]

    return xr.concat(flux_stats, dim="time")


def make_concentration_outputs(inv_outs: list[InversionOutput]) -> xr.Dataset:
    predictive_vars = ["y", "mu_bc"]
    preds_list = [inv_out.get_trace_dataset(var_names=predictive_vars) for inv_out in inv_outs]

    conc_stats = []
    for preds in preds_list:
        var_names = ["mu_bc_posterior", "mu_bc_prior"]
        stats = calculate_stats(preds, name="Y", chunk_dim="nmeasure", var_names=var_names, report_mode=True)

        var_names = ["y_posterior_predictive", "y_prior_predictive"]
        stats.extend(
            calculate_stats(
                preds,
                name="Y",
                chunk_dim="nmeasure",
                var_names=var_names,
                report_mode=True,
                add_bc_suffix=True,
            )
        )

        conc_stats.append(xr.merge(stats).unstack("nmeasure"))

    conc_output = xr.concat(conc_stats, dim="time")

    return conc_output


def main(
    species: str,
    output_file_path: str,
    country_files_root: str,
    min_model_error: float = 0.1,
    n_files: Optional[int] = None,
) -> tuple[xr.Dataset, xr.Dataset]:
    inv_outs = get_inversion_outputs_with_samples(
        species=species, output_file_path=output_file_path, min_model_error=min_model_error
    )

    if n_files is not None:
        inv_outs = inv_outs[:n_files]

    # make country and flux output
    country_output = make_country_output(inv_outs, country_files_root)
    flux_output = make_flux_outputs(inv_outs)

    template_file = str(PARIS_FORMATTING_PATH / "PARIS_Lagrangian_inversion_flux_EUROPE.cdl")
    emissions_attrs = get_data_var_attrs(template_file, species)

    # renaming as in latest .cdl file from Stephan
    rename_dict = {"lat": "latitude",
                   "lon": "longitude",
                   "probs": "percentile"}

    vars_to_drop = []

    for dv in country_output.data_vars:
        if str(dv).startswith("q"):
            if str(dv).endswith("apost"):
                rename_dict[str(dv)] = "percentile_country_flux_total_posterior"
            else:
                rename_dict[str(dv)] = "percentile_country_flux_total_prior"
        elif str(dv).endswith("_mode"):
            if str(dv)[:-5].endswith("apost"):
                rename_dict[str(dv)] = "country_flux_total_posterior"
            else:
                rename_dict[str(dv)] = "country_flux_total_prior"
        else:
            vars_to_drop.append(dv)  # drop means

    for dv in flux_output.data_vars:
        if str(dv).startswith("q"):
            if str(dv).endswith("apost"):
                rename_dict[str(dv)] = "percentile_flux_total_posterior"
            else:
                rename_dict[str(dv)] = "percentile_flux_total_prior"
        elif str(dv).endswith("_mode"):
            if str(dv)[:-5].endswith("apost"):
                rename_dict[str(dv)] = "flux_total_posterior"
            else:
                rename_dict[str(dv)] = "flux_total_prior"
        else:
            vars_to_drop.append(dv)  # drop means

    # merge and process names, attrs
    emissions = (
        xr.merge([flux_output, country_output])
        .pipe(convert_time_to_unix_epoch)
        .drop_vars(vars_to_drop)
        .rename(rename_dict)
        .pipe(add_variable_attrs, emissions_attrs)
    )

    emissions.attrs = make_global_attrs("flux")


    # make concentration outputs
    conc_output = make_concentration_outputs(inv_outs)
    y_obs = xr.concat([inv_out.get_obs().unstack("nmeasure") for inv_out in inv_outs], dim="time")

    conc_attrs = get_data_var_attrs(
        str(PARIS_FORMATTING_PATH / "netcdf_template_concentrations_bm_edits.txt")
    )

    units = float(y_obs.attrs["units"].split(" ")[0])  # e.g. get 1e-12 from "1e-12 mol/mol"

    concentrations = (
        xr.merge([y_obs, conc_output])
        .pipe(convert_time_to_unix_epoch)
        .rename(probs="quantile")
        .pipe(add_variable_attrs, conc_attrs, units)
    )

    concentrations.attrs = make_global_attrs("conc")

    return emissions, concentrations
