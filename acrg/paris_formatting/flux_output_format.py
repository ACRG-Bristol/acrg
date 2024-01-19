#!/usr/bin/env python
import json
from typing import Any, Optional

import xarray as xr
import pandas as pd
import numpy as np
import pymc as pm
import arviz as az

from attribute_parsers import make_global_attrs


def parseprior(name: str, prior_params: dict[str, Any], **kwargs):
    """
    Parses all continuous distributions for PyMC 3.8:
    https://docs.pymc.io/api/distributions/continuous.html
    This format requires updating when the PyMC distributions update,
    but is safest for code execution
    -----------------------------------
    Args:
      name (str):
        name of variable in the pymc model
      prior_params (dict):
        dict of parameters for the distribution,
        including 'pdf' for the distribution to use
      shape (array):
        shape of distribution to be created.
        Default shape = () is the same as used by PyMC3
    -----------------------------------
    """
    functiondict = {
        "uniform": pm.Uniform,
        "flat": pm.Flat,
        "halfflat": pm.HalfFlat,
        "normal": pm.Normal,
        "truncatednormal": pm.TruncatedNormal,
        "halfnormal": pm.HalfNormal,
        "skewnormal": pm.SkewNormal,
        "beta": pm.Beta,
        "kumaraswamy": pm.Kumaraswamy,
        "exponential": pm.Exponential,
        "laplace": pm.Laplace,
        "studentt": pm.StudentT,
        "halfstudentt": pm.HalfStudentT,
        "cauchy": pm.Cauchy,
        "halfcauchy": pm.HalfCauchy,
        "gamma": pm.Gamma,
        "inversegamma": pm.InverseGamma,
        "weibull": pm.Weibull,
        "lognormal": pm.Lognormal,
        "chisquared": pm.ChiSquared,
        "wald": pm.Wald,
        "pareto": pm.Pareto,
        "exgaussian": pm.ExGaussian,
        "vonmises": pm.VonMises,
        "triangular": pm.Triangular,
        "gumbel": pm.Gumbel,
        "rice": pm.Rice,
        "logistic": pm.Logistic,
        "logitnormal": pm.LogitNormal,
        "interpolated": pm.Interpolated,
    }

    pdf = prior_params.pop("pdf")
    return functiondict[pdf.lower()](name, **prior_params, **kwargs)


def get_prior_samples(
    ds: xr.Dataset,
    xprior: Optional[dict[str, Any]] = None,
    bcprior: Optional[dict[str, Any]] = None,
    sigprior: Optional[dict[str, Any]] = None,
    min_model_error: float = 20.0,
    coord_dim: str = "nmeasure",
) -> az.InferenceData:
    """Sample prior and prior predictive distributions from RHIME model."""
    if xprior is None:
        xprior = {"pdf": "truncatednormal", "mu": 1.0, "sigma": 1.0, "lower": 0.0}
    if bcprior is None:
        bcprior = {"pdf": "truncatednormal", "mu": 1.0, "sigma": 0.1, "lower": 0.0}
    if sigprior is None:
        sigprior = {"pdf": "uniform", "lower": 0.0, "upper": 1.0}

    coords_dict = {coord_dim: ds[coord_dim]}
    for dim in ["nlatent", "nBC", "nsigma_site", "nsigma_time"]:
        coords_dict[dim] = ds[dim]

    with pm.Model(coords=coords_dict) as model:
        x = parseprior("x", xprior, dims="nlatent")
        bc = parseprior("bc", bcprior, dims="nBC")
        sigma = parseprior("sigma", sigprior, dims=("nsigma_site", "nsigma_time"))

        muBC = pm.math.dot(ds.bcsensitivity.values, bc)
        mu = pm.math.dot(ds.xsensitivity.values, x) + muBC

        mult_error = (
            np.abs(pm.math.dot(ds.xsensitivity.values, x))
            * sigma[ds.siteindicator.values.astype(int), ds.sigmafreqindex.values.astype(int)]
        )
        epsilon = pm.math.sqrt(ds.Yerror.values**2 + mult_error**2 + min_model_error**2)
        ymodbc = pm.Normal("ymodbc", muBC, epsilon, dims=coord_dim)
        ymod = pm.Normal("ymod", mu, epsilon, dims=coord_dim)
        idata = pm.sample_prior_predictive(1000)

    return idata


def format_concentrations(ds0: xr.Dataset) -> xr.Dataset:
    """Return the formatted concentration time series given hbmcmc outputs.

    TODO: split out quantile function
    TODO: move sampling to inversions
    """
    ## rename coordinates
    ds0 = ds0.rename_vars(
        stepnum="steps", paramnum="nlatent", numBC="nBC", measurenum="nmeasure", UInum="nUI"
    )

    # Concentrations (timeseries) netCDF

    # sample prior predictive
    min_model_error = 20.0
    idata = get_prior_samples(ds0)

    # make quantile vars
    probs = [0.025, 0.159, 0.841, 0.975]

    # qYmod
    y_mod = idata.prior.ymod.isel(chain=0)
    y_mod_q = np.quantile(y_mod, probs, axis=0)

    qYmod_da = xr.DataArray(y_mod_q.T, coords=[ds0.nmeasure, probs], dims=["nmeasure", "probs"], name="qYmod")

    # uYaprioriBC
    y_modbc = idata.prior.ymodbc.isel(chain=0)
    y_aprioribc_u = 0.5 * np.diff(np.quantile(y_modbc, probs[1:3], axis=0), axis=0).squeeze()
    ds0["uYaprioriBC"] = xr.DataArray(
        y_aprioribc_u, coords=[ds0.nmeasure], dims=["nmeasure"], name="qYaprioriBC"
    )

    # We need to combine the 68 and 95 confidence intervals into a single variable.
    y_apost_q = np.vstack(
        [ds0.Ymod95.values[:, 0], ds0.Ymod68.values[:, 0], ds0.Ymod68.values[:, 1], ds0.Ymod95.values[:, 1]]
    )
    qYapost_da = xr.DataArray(
        y_apost_q.T, coords=[ds0.nmeasure, probs], dims=["nmeasure", "probs"], name="qYapost"
    )

    y_apostbc_q = np.vstack(
        [
            ds0.Ymod95BC.values[:, 0],
            ds0.Ymod68BC.values[:, 0],
            ds0.Ymod68BC.values[:, 1],
            ds0.Ymod95BC.values[:, 1],
        ]
    )
    qYapostBC_da = xr.DataArray(
        y_apostbc_q.T, coords=[ds0.nmeasure, probs], dims=["nmeasure", "probs"], name="qYapostBC"
    )

    ds0 = ds0.drop_vars(["Ymod68", "Ymod95", "Ymod68BC", "Ymod95BC"])
    ds0 = ds0.drop_dims("nUI")

    ds0 = xr.merge([ds0, qYmod_da, qYapost_da, qYapostBC_da])

    # rename other timeseries variables
    ds0 = ds0.rename_vars(Yerror="uYobs", Ymodmean="Yapost", YmodmeanBC="YapostBC")

    # Select timeseries data
    time_vars = [dv for dv in ds0.data_vars if "nmeasure" in list(ds0[dv].coords)]
    time_vars = [
        tv
        for tv in time_vars
        if all(x not in str(tv) for x in ["mode", "median", "off", "sensitivity", "sigma"])
    ]

    ds = ds0[time_vars]

    #  `nmeasure` to `(time, nsite)`
    site_nums = np.arange(ds.siteindicator.max() + 1)

    ds = (
        xr.concat(
            [
                ds.where(ds.siteindicator == site_num, drop=True)
                .expand_dims({"nsite": [site_num]})
                .assign_coords(nmeasure=ds.Ytime.where(ds.siteindicator == site_num, drop=True))
                .rename_vars(nmeasure="time")
                for site_num in site_nums
            ],
            dim="nsite",
        )
        # .swap_dims(nmeasure="time")
        .drop_vars(["Ytime", "siteindicator"])
    )

    # extend lon, lat along time dimension
    def nsites_to_nsite_nmeasure(da: xr.DataArray) -> xr.DataArray:
        da = xr.DataArray(da.values, coords=[ds.nsite], dims=["nsite"])
        da, _ = xr.broadcast(da, ds.Yobs)
        return da

    ds = ds.merge(
        {
            "sitenames": nsites_to_nsite_nmeasure(ds0.sitenames),
            "longitude": nsites_to_nsite_nmeasure(ds0.sitelons),
            "latitude": nsites_to_nsite_nmeasure(ds0.sitelats),
        }
    )
    ds = ds.swap_dims(nmeasure="time")

    # Updating attributes
    # Attributes for data vars, extracted from netcdf template for concentrations:

    with open("data_var_attrs.json", "r") as f:
        data_var_attrs = json.load(f)

    # Some attributes will only be correct if we change the units of the values in the data variable:
    # - `time` needs to be changed to seconds from UNIX epoch
    # - all `Y` varaibles need to be multiplied by `1e-9` to be in mol/mol instead of ppb
    #
    # The fill value variables are probably set by `to_netcdf`?
    data_var_attrs = {
        k: {j: w for j, w in v.items() if "FillValue" not in j} for k, v in data_var_attrs.items()
    }

    # Formatting `time` coord to UNIX epoch
    ds = ds.assign_coords(time=(pd.DatetimeIndex(ds.time) - pd.Timestamp("1970-01-01")) // pd.Timedelta("1s"))

    # Adding attributes
    current_units = 1e-9  # mol/mol
    for k, v in data_var_attrs.items():
        if k in ds.data_vars:
            if "units" in v and "mol/mol" in v["units"]:
                ds[k] = current_units * ds[k]
            ds[k].attrs = v
        elif k in ds.coords:
            ds.coords[k].attrs = v

    # Updating global attributes
    ds.attrs = make_global_attrs("conc")

    return ds
