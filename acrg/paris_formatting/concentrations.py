from typing import Union

import numpy as np
import xarray as xr

from attribute_parsers import (
    add_variable_attrs,
    conc_template_path,
    convert_time_to_unix_epoch,
    get_data_var_attrs,
    make_global_attrs,
)
from process_rhime_output import InversionOutput
from stats import calculate_stats


def _make_concentration_outputs(inv_outs: list[InversionOutput], report_mode: bool = False) -> xr.Dataset:
    """Extract y and mu_bc traces, calculate stats, and combine along time axis."""
    use_bc = "mu_bc" in inv_outs[0].trace["prior"].data_vars

    predictive_vars = ["y", "mu_bc"] if use_bc else ["y"]
    preds_list = [inv_out.get_trace_dataset(var_names=predictive_vars) for inv_out in inv_outs]

    conc_stats = []
    for preds in preds_list:
        if use_bc:
            var_names = ["mu_bc_posterior", "mu_bc_prior"]
            stats = calculate_stats(
                preds,
                name="Y",
                chunk_dim="nmeasure",
                chunk_size=1,
                var_names=var_names,
                report_mode=report_mode,
                add_bc_suffix=True,
            )
        else:
            stats = []

        var_names = ["y_posterior_predictive", "y_prior_predictive"]
        stats.extend(
            calculate_stats(
                preds,
                name="Y",
                chunk_dim="nmeasure",
                var_names=var_names,
                report_mode=report_mode,
            )
        )

        conc_stats.append(xr.merge(stats).unstack("nmeasure"))

    conc_output = xr.concat(conc_stats, dim="time")

    return conc_output


def shift_measurement_time_to_midpoint(
    ds: xr.Dataset, period: str = "4h"
) -> xr.Dataset:
    """Adjust `time` coordinate of concentrations to represent half averaging "period"."""
    time_shifted = ds["time"].astype("datetime64[ns]") + np.timedelta64(int(period[:-1]), period[-1]) / 2
    ds = ds.assign_coords(time=time_shifted)
    return ds


def make_concentration_outputs(
    inv_outs: list[InversionOutput], report_mode: bool = False, avr_obs_period: str = "4h"
) -> xr.Dataset:
    conc_output = _make_concentration_outputs(inv_outs, report_mode=report_mode)
    y_obs = xr.concat([inv_out.get_obs() for inv_out in inv_outs], dim="time")
    y_obs_err = xr.concat([inv_out.get_obs_err() for inv_out in inv_outs], dim="time")

    conc_attrs = get_data_var_attrs(conc_template_path)

    units = float(y_obs.attrs["units"].split(" ")[0])  # e.g. get 1e-12 from "1e-12 mol/mol"

    # renaming as in latest .cdl file from Stephan
    rename_dict_conc = {"probs": "percentile", "site": "nsite", "Yerror": "uYobs", "qYapriori": "qYmod"}
    if "qYapostBC" in conc_output.data_vars:
        vars_to_drop_conc = ["qYapostBC", "qYaprioriBC"]
    else:
        vars_to_drop_conc = []

    # merge and process names, attrs
    concentrations = (
        xr.merge([y_obs, y_obs_err, conc_output])
        .pipe(shift_measurement_time_to_midpoint, avr_obs_period)
        .pipe(convert_time_to_unix_epoch, "1s")
        .drop_vars(vars_to_drop_conc)
        .rename(rename_dict_conc)
        .pipe(add_variable_attrs, conc_attrs, units)
        .transpose("time", "percentile", "nsite")
    )

    # add sitenames variable and remove nsite coordinate
    # TODO: should sitenames be a coordinate?
    # it would be better to leave site as a dimension coordinate so we can just do `.sel(site="MHD")` etc.
    concentrations = concentrations.assign(sitenames=("nsite", concentrations.nsite.values.astype("|S3")))
    concentrations["sitenames"].attrs["long_name"] = "identifier of site"
    del concentrations["nsite"]

    concentrations.attrs = make_global_attrs("conc")

    return concentrations
