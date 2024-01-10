#!/usr/bin/env python
import json
import xarray as xr
import pandas as pd
import numpy as np
import pymc as pm
import arviz as az


from attribute_parsers import make_global_attrs


def get_prior_samples(
    ds: xr.Dataset, min_model_error: float = 20.0, coord_dim: str = "nmeasure"
) -> az.InferenceData:
    with pm.Model(coords={coord_dim: ds[coord_dim]}) as model:
        x = pm.TruncatedNormal("x", mu=1.0, sigma=1.0, lower=0.0, shape=ds.dims["nlatent"])
        bc = pm.TruncatedNormal("bc", mu=1.0, sigma=0.1, lower=0.0, shape=ds.dims["nBC"])
        sigma = pm.Uniform(
            "sigma", lower=0.1, upper=1, shape=(ds.dims["nsigma_site"], ds.dims["nsigma_time"])
        )

        muBC = pm.Deterministic("muBC", pm.math.dot(ds.bcsensitivity.values, bc), dims=coord_dim)
        mu = pm.Deterministic("mu", pm.math.dot(ds.xsensitivity.values, x) + muBC, dims=coord_dim)

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
