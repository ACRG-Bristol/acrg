"""
Script for creating PARIS outputs from RHIME inversion outputs.
"""
from argparse import ArgumentParser
from functools import partial
from pathlib import Path
from typing import Literal, Optional, Union

import numpy as np
import xarray as xr
from openghg.util import timestamp_now

from attribute_parsers import (
    add_variable_attrs,
    get_data_var_attrs,
    make_global_attrs,
    get_country_code,
    convert_time_to_unix_epoch,
)
from countries import Countries
from array_ops import sparse_xr_dot
from process_rhime_output import InversionOutput
from stats import calculate_stats


paris_formatting_path = Path(__file__).parent


def get_netcdf_files(directory: Union[str, Path], filename_search: Optional[str] = None) -> list[Path]:
    """Get list of paths to netCDF files, optionally filtered by `filename_search` string."""
    if isinstance(directory, str):
        directory = Path(directory)

    if filename_search is None:
        file_list = sorted(directory.glob("*.nc"))
    else:
        file_list = sorted(directory.glob(f"*{filename_search}*.nc"))

    return file_list


def get_inversion_outputs_with_samples(
    species: str, output_file_path: str, min_model_error: float = 0.1
) -> list[InversionOutput]:
    """Create a list of InversionOutputs given a path to RHIME inversion outputs."""
    files = get_netcdf_files(output_file_path, filename_search=species.upper())
    inv_outs = [InversionOutput.from_rhime(xr.open_dataset(file), min_model_error) for file in files]

    for inv_out in inv_outs:
        inv_out.sample_predictive_distributions()

    return inv_outs


def get_time_point(inv_out: InversionOutput, time_point: Literal["start", "midpoint"]) -> np.datetime64:
    """Get time point to represent inversion period."""
    if time_point == "start":
        return inv_out.start_time()

    if time_point == "midpoint":
        return inv_out.period_midpoint()

    raise ValueError(f"time_point must be 'start' or 'midpoint'; given {time_point}.")


def make_country_output(
        inv_outs: list[InversionOutput], country_files_root: str, code: Literal["alpha2", "alpha3"] = "alpha3", time_point: Literal["start", "midpoint"] = "midpoint"
) -> xr.Dataset:
    """Calculate country stats."""
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

    time_func = partial(get_time_point, time_point=time_point)

    country_traces = [
        countries.get_country_trace("sf6", inv_out).expand_dims({"time": [time_func(inv_out)]})
        for inv_out in inv_outs
    ]

    country_traces_merged = xr.concat(country_traces, dim="time")

    country_output = xr.merge(
        calculate_stats(country_traces_merged, "country", chunk_dim="country", chunk_size=1)
    )

    # TODO: calculate stats for composite regions e.g. BE-NE-LX
    # to do this, we need country codes applied to country_traces_merged

    # apply `get_country_code` to each element of `country` coordinate
    country_codes = list(map(partial(get_country_code, code=code), map(str, country_output.country.values)))
    country_output = country_output.assign_coords(country=country_codes)

    return country_output


def make_flux_outputs(inv_outs: list[InversionOutput], time_point: Literal["start", "midpoint"] = "midpoint") -> xr.Dataset:
    """Make flux output dataset"""

    # calculate stats on flux traces
    traces = [inv_out.get_trace_dataset(convert_nmeasure=False, var_names="x") for inv_out in inv_outs]
    stats = [xr.merge(calculate_stats(trace, "flux", chunk_dim="nx")) for trace in traces]

    time_func = partial(get_time_point, time_point=time_point)

    # multiply stats by matrix mapping basis regions to lat/lon
    flux_stats = [
        sparse_xr_dot((inv_out.flux * inv_out.basis), stats_ds).expand_dims({"time": [time_func(inv_out)]})
        for inv_out, stats_ds in zip(inv_outs, stats)
    ]

    return xr.concat(flux_stats, dim="time")


def make_concentration_outputs(inv_outs: list[InversionOutput]) -> xr.Dataset:
    """Extract y and mu_bc traces, calculate stats, and combine along time axis."""
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


def rename_drop_dvs_for_template(ds: xr.Dataset, var_name: str) -> tuple[dict[str, str], list[str]]:
    """Returns dict to rename data vars and list of data vars to drop."""
    rename_dict = {}
    vars_to_drop = []

    for dv in ds.data_vars:
        if str(dv).startswith("q"):
            if str(dv).endswith("apost"):
                rename_dict[str(dv)] = f"percentile_{var_name}_total_posterior"
            else:
                rename_dict[str(dv)] = f"percentile_{var_name}_total_prior"
        elif str(dv).endswith("_mode"):
            if str(dv)[:-5].endswith("apost"):
                rename_dict[str(dv)] = f"{var_name}_total_posterior"
            else:
                rename_dict[str(dv)] = f"{var_name}_total_prior"
        else:
            vars_to_drop.append(dv)  # drop means

    return rename_dict, vars_to_drop


def main(
    species: str,
    output_file_path: str,
    country_files_root: str,
    min_model_error: float,
    n_files: Optional[int] = None,
) -> tuple[xr.Dataset, xr.Dataset]:
    """Create formatted PARIS emissions and concentrations datasets.

    Args:
        species: species used in inversion
        output_file_path: path to directory containing RHIME outputs
        country_files_root: path to directory containing country files
        min_model_error: the minimum model error used with the inversion
        n_files: number of output files to process. This is mainly to keep runs small for testing.

    Returns:
        emissions dataset and concentrations dataset
    """
    inv_outs = get_inversion_outputs_with_samples(
        species=species, output_file_path=output_file_path, min_model_error=min_model_error
    )

    if n_files is not None:
        inv_outs = inv_outs[:n_files]

    # make country and flux output
    country_output = make_country_output(inv_outs, country_files_root)
    flux_output = make_flux_outputs(inv_outs)

    template_file = str(paris_formatting_path / "PARIS_Lagrangian_inversion_flux_EUROPE.cdl")
    emissions_attrs = get_data_var_attrs(template_file, species)

    # renaming as in latest .cdl file from Stephan
    rename_dict = {"lat": "latitude", "lon": "longitude", "probs": "percentile"}

    vars_to_drop = []

    rename_dict_country, vars_to_drop_country = rename_drop_dvs_for_template(country_output, "country_flux")
    rename_dict_flux, vars_to_drop_flux = rename_drop_dvs_for_template(flux_output, "flux")

    rename_dict.update(rename_dict_country)
    rename_dict.update(rename_dict_flux)

    vars_to_drop.append(vars_to_drop_country)
    vars_to_drop.append(vars_to_drop_flux)

    # merge and process names, attrs
    emissions = (
        xr.merge([flux_output, country_output])
        .pipe(convert_time_to_unix_epoch, "1D")
        .drop_vars(vars_to_drop)
        .rename(rename_dict)
        .pipe(add_variable_attrs, emissions_attrs)
    )

    emissions.attrs = make_global_attrs("flux")

    # make concentration outputs
    conc_output = make_concentration_outputs(inv_outs)
    y_obs = xr.concat([inv_out.get_obs() for inv_out in inv_outs], dim="time")

    conc_attrs = get_data_var_attrs(
        str(paris_formatting_path / "netcdf_template_concentrations_bm_edits.txt")
    )

    units = float(y_obs.attrs["units"].split(" ")[0])  # e.g. get 1e-12 from "1e-12 mol/mol"

    concentrations = (
        xr.merge([y_obs, conc_output])
        .pipe(convert_time_to_unix_epoch, "1D")
        .rename(probs="quantile")
        .pipe(add_variable_attrs, conc_attrs, units)
    )

    concentrations.attrs = make_global_attrs("conc")

    return emissions, concentrations


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-s", "--species", type=str, help="species used in inversion; will be used to filter files"
    )
    parser.add_argument("-r", "--rhime-outputs-path", type=str, help="path to RHIME outputs")
    parser.add_argument(
        "-c",
        "--country-files-path",
        type=str,
        help="path to country files; must contain 'country_EUROPE.nc' and 'country-ukmo_EUROPE.nc'",
    )
    parser.add_argument("-m", "--min-model-error", type=float, help="min. model error used in inversion")
    parser.add_argument("-o", "--output-path", type=str, help="path to dir to write formatted outputs")
    parser.add_argument("-t", "--output-tag", type=str, help="tag to add to output file names")

    args = parser.parse_args()

    emissions, concentrations = main(
        species=args.species,
        output_file_path=args.rhime_outputs_path,
        country_files_root=args.country_files_path,
        min_model_error=args.min_model_error,
    )

    if args.tag:
        tag = args.tag
    else:
        date, time = str(timestamp_now()).split(" ")
        tag = date + "_" + time.split(".")[0]

    output_path = Path(args.output_path)
    conc_output_path = output_path / f"PARIS_concentrations_{args.species}_{tag}.nc"
    emissions_output_path = output_path / f"PARIS_emissions_{args.species}_{tag}.nc"
