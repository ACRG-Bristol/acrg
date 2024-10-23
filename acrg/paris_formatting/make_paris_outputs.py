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
    convert_time_to_unix_epoch,
    flux_template_path,
    get_country_code,
    get_data_var_attrs,
    make_global_attrs,
)
from concentrations import make_concentration_outputs
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
    species: str,
    output_file_path: str,
    ndraw: int = 10000,
    n_files: Optional[int] = None,
    pol_from_obs: bool = False,
    no_model_error: bool = False,
) -> list[InversionOutput]:
    """Create a list of InversionOutputs given a path to RHIME inversion outputs."""
    files = get_netcdf_files(output_file_path, filename_search=species.upper())

    if n_files is not None:
        files = files[:n_files]

    inv_outs = [
        InversionOutput.from_rhime(
            xr.open_dataset(file), pol_from_obs=pol_from_obs, no_model_error=no_model_error, ndraw=ndraw
        )
        for file in files
    ]

    for inv_out in inv_outs:
        inv_out.sample_predictive_distributions(ndraw=ndraw)

    return inv_outs


def get_time_point(inv_out: InversionOutput, time_point: Literal["start", "midpoint"]) -> np.datetime64:
    """Get time point to represent inversion period."""
    if time_point == "start":
        return inv_out.start_time()

    if time_point == "midpoint":
        return inv_out.period_midpoint()

    raise ValueError(f"time_point must be 'start' or 'midpoint'; given {time_point}.")


def make_country_output(
    species: str,
    inv_outs: list[InversionOutput],
    countries: Countries,
    code: Literal["alpha2", "alpha3"] = "alpha3",
    time_point: Literal["start", "midpoint"] = "midpoint",
    report_mode: bool = True,
) -> xr.Dataset:
    """Calculate country stats."""
    time_func = partial(get_time_point, time_point=time_point)

    country_traces = [
        countries.get_country_trace(species, inv_out).expand_dims({"time": [time_func(inv_out)]})
        for inv_out in inv_outs
    ]

    country_traces_merged = xr.concat(country_traces, dim="time")

    # apply `get_country_code` to each element of `country` coordinate
    country_codes = list(
        map(partial(get_country_code, code=code), map(str, country_traces_merged.country.values))
    )
    country_traces_merged = country_traces_merged.assign_coords(country=country_codes)

    # add country regions
    regions_dict = {
        "BELUX": "BEL-LUX",
        "BENELUX": "BEL-LUX-NLD",
        "CW_EU": "AUT-BEL-CHE-CZE-DEU-ESP-FRA-GBR-HRV-HUN-IRL-ITA-LUX-NLD-POL-PRT-SVK-SVN",
        "EU_GRP2": "AUT-BEL-CHE-DEU-DNK-FRA-GBR-IRL-ITA-LUX-NLD",
        "NW_EU": "BEL-DEU-DNK-FRA-GBR-IRL-LUX-NLD",
        "NW_EU2": "BEL-DEU-FRA-GBR-IRL-LUX-NLD",
        "NW_EU_CONTINENT": "BEL-DEU-FRA-LUX-NLD",
    }

    try:
        region_traces = []
        for region, countries_str in regions_dict.items():
            region_countries = countries_str.split("-")
            region_ds = (
                country_traces_merged.sel(country=region_countries).sum("country").expand_dims({"country": [region]})
            )
            region_traces.append(region_ds)
    except KeyError as e:
        print(f"A KeyError occurred: {e}")
        print(f"Warning: Country regions ({', '.join(regions_dict.keys())}) were not added.")
        all_traces_merged = country_traces_merged
    else:
        region_traces_merged = xr.concat(region_traces, dim="country")
        # combine country and region traces
        all_traces_merged = xr.merge([country_traces_merged, region_traces_merged])

    country_output = xr.merge(
        calculate_stats(
            all_traces_merged, "country", chunk_dim="country", chunk_size=1, report_mode=report_mode
        )
    )

    return country_output


def make_flux_outputs(
    inv_outs: list[InversionOutput],
    time_point: Literal["start", "midpoint"] = "midpoint",
    report_mode: bool = False,
    inversion_grid: bool = False,
) -> xr.Dataset:
    """Make flux output dataset"""

    # calculate stats on flux traces
    traces = [inv_out.get_trace_dataset(convert_nmeasure=False, var_names="x") for inv_out in inv_outs]
    stats = [
        xr.merge(calculate_stats(trace, "flux", chunk_dim="nx", report_mode=report_mode)) for trace in traces
    ]

    time_func = partial(get_time_point, time_point=time_point)

    if inversion_grid is False:
        # multiply stats by matrix mapping basis regions to lat/lon
        flux_stats = [
            sparse_xr_dot((inv_out.flux * inv_out.basis), stats_ds).expand_dims({"time": [time_func(inv_out)]})
            for inv_out, stats_ds in zip(inv_outs, stats)
        ]
    else:
        # sum prior flux over basis regions
        agg_fluxes = [sparse_xr_dot(inv_out.basis, inv_out.flux) for inv_out in inv_outs]

        flux_stats = [
            sparse_xr_dot(inv_out.basis, agg_flux * stats_ds).expand_dims({"time": [time_func(inv_out)]})
            for inv_out, agg_flux, stats_ds in zip(inv_outs, agg_fluxes, stats)
        ]

    return xr.concat(flux_stats, dim="time")


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
        else:
            if str(dv).endswith("apost"):
                rename_dict[str(dv)] = f"{var_name}_total_posterior"
            else:
                rename_dict[str(dv)] = f"{var_name}_total_prior"

    return rename_dict, vars_to_drop


def main(
    species: str,
    output_file_path: str,
    country_file_path: str,
    avr_obs_period: str,
    n_files: Optional[int] = None,
    return_concentrations: bool = True,
    report_mf_mode: bool = False,
    report_em_mode: bool = False,
    pol_from_obs: bool = False,
    no_model_error: bool = False,
    ndraw: int = 10000,
) -> tuple[xr.Dataset, Optional[xr.Dataset]]:
    """Create formatted PARIS emissions and concentrations datasets.

    Args:
        species: species used in inversion
        output_file_path: path to directory containing RHIME outputs
        country_files_root: path to directory containing country files
        avr_obs_period: the averaging period for measurements used in inversion
        n_files: number of output files to process. This is mainly to keep runs small for testing.
        return_concentrations: if False, only country and flux outputs are returned. (None is returned for concentrations.)
        report_mf_mode: if True, use mode for concentration prior/posterior predictives (i.e. for y and y BC)
        report_em_mode: if True, use mode for country and flux prior/posterior totals
        pol_from_obs: if True, this means that pollution_events_from_obs=True was specified

    Returns:
        emissions dataset and concentrations dataset
    """
    inv_outs = get_inversion_outputs_with_samples(
        species=species,
        output_file_path=output_file_path,
        n_files=n_files,
        pol_from_obs=pol_from_obs,
        no_model_error=no_model_error,
        ndraw=ndraw,
    )

    # make country and flux output
    countries = Countries(xr.open_dataset(country_file_path))

    country_output = make_country_output(species, inv_outs, countries, report_mode=report_em_mode)
    flux_output = make_flux_outputs(inv_outs, report_mode=report_em_mode)

    emissions_attrs = get_data_var_attrs(flux_template_path, species)

    # renaming as in latest .cdl file from Stephan
    rename_dict = {"lat": "latitude", "lon": "longitude", "probs": "percentile"}

    vars_to_drop = []

    rename_dict_country, vars_to_drop_country = rename_drop_dvs_for_template(country_output, "country_flux")
    rename_dict_flux, vars_to_drop_flux = rename_drop_dvs_for_template(flux_output, "flux")

    rename_dict.update(rename_dict_country)
    rename_dict.update(rename_dict_flux)

    vars_to_drop.extend(vars_to_drop_country)
    vars_to_drop.extend(vars_to_drop_flux)

    # merge and process names, attrs
    emissions = (
        xr.merge([flux_output, country_output * 1e-3])  # convert g/yr to kg/yr
        .pipe(convert_time_to_unix_epoch, "1s")
        .drop_vars(vars_to_drop)
        .rename(rename_dict)
        .pipe(add_variable_attrs, emissions_attrs)
        .transpose("time", "latitude", "longitude", "percentile", "country")
    )

    # add flux on inversion grid
    flux_output_inversion_grid = make_flux_outputs(inv_outs, report_mode=report_em_mode, inversion_grid=True)
    rename_dict_flux_inversion_grid, vars_to_drop_flux_inversion_grid = rename_drop_dvs_for_template(flux_output_inversion_grid, "flux")
    rename_dict_flux_inversion_grid2 = {v: f"{v}_inversion_grid" for v in rename_dict_flux_inversion_grid.values()}

    emissions_inversion_grid = (
        flux_output_inversion_grid
        .pipe(convert_time_to_unix_epoch, "1s")
        .drop_vars(vars_to_drop_flux_inversion_grid)
        .rename(rename_dict_flux_inversion_grid)
        .pipe(add_variable_attrs, emissions_attrs)
        .transpose("time", "latitude", "longitude", "percentile", "country")
        .rename(rename_dict_flux_inversion_grid2)  # update name to have _inversion_grid at the end; do this here to get attributes from old template
    )

    emissions = emissions.merge(emissions_inversion_grid)

    emissions.attrs = make_global_attrs("flux")

    # make concentration outputs
    if return_concentrations is True:
        concentrations = make_concentration_outputs(inv_outs, report_mode=report_mf_mode, avr_obs_period=avr_obs_period)
    else:
        concentrations = None

    return emissions, concentrations


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "-s", "--species", type=str, help="species used in inversion; will be used to filter files"
    )
    parser.add_argument("-r", "--rhime-outputs-path", type=str, help="path to RHIME outputs")
    parser.add_argument(
        "-c",
        "--country-file-path",
        type=str,
        help="path to country file; e.g. '/group/acrg/chem/LPDM/countries/country_EUROPE.nc'",
    )
    parser.add_argument(
        "-p", "--avr-obs-period", type=str, help="averaging period for measurements used in inversion"
    )
    parser.add_argument("-o", "--output-path", type=str, help="path to dir to write formatted outputs")
    parser.add_argument("-t", "--output-tag", type=str, help="tag to add to output file names")
    parser.add_argument("-n", "--n-files", type=int, help="number of files to process")
    parser.add_argument(
        "--no-conc",
        action="store_true",
        default=False,
        help="if set, only process emissions and country totals.",
    )
    parser.add_argument(
        "--mode",
        action="store_true",
        default=False,
        help="if set, report mode for concentrations/mole fractions (by default, mean is reported).",
    )
    parser.add_argument(
        "--em-mode",
        action="store_true",
        default=False,
        help="if set, report mode for country and flux totals (by default, mean is reported).",
    )
    parser.add_argument(
        "--pol-obs",
        action="store_true",
        default=False,
        help="if set, this means that pollution_events_from_obs=True was specified",
    )
    parser.add_argument(
        "--no-model-error",
        action="store_true",
        default=False,
        help="if set, this means that no_model_error=True was specified",
    )
    parser.add_argument("--ndraw", type=int, default=10000, help="number of prior/predictive samples to produce")

    args = parser.parse_args()

    emissions, concentrations = main(
        species=args.species,
        output_file_path=args.rhime_outputs_path,
        country_file_path=args.country_file_path,
        avr_obs_period=args.avr_obs_period,
        n_files=args.n_files,
        return_concentrations=(not args.no_conc),
        report_mf_mode=args.mode,
        report_em_mode=args.em_mode,
        pol_from_obs=args.pol_obs,
        no_model_error=args.no_model_error,
        ndraw=args.ndraw,
    )

    output_path = Path(args.output_path)

    if not output_path.exists():
        output_path.mkdir(parents=True)

    if args.output_tag:
        tag = args.output_tag
        emissions_output_path = output_path / f"{tag}.nc"
        conc_output_path = output_path / f"{tag}_concentrations.nc"
    else:
        date, time = str(timestamp_now()).split(" ")
        tag = date + "_" + time.split(".")[0].replace(":", "")
        emissions_output_path = output_path / f"PARIS_emissions_{args.species}_{tag}.nc"
        conc_output_path = output_path / f"PARIS_concentrations_{args.species}_{tag}.nc"

    emissions.to_netcdf(emissions_output_path)

    if not args.no_conc:
        concentrations.to_netcdf(conc_output_path)  # type: ignore
