"""
Code to process RHIME outputs into format that is easier to
manipulate.
"""
from dataclasses import dataclass
from typing import Optional, TypeVar, Union

import arviz as az
import numpy as np
import pymc as pm
import xarray as xr

from array_ops import get_xr_dummies
from sampling import (
    convert_idata_to_dataset,
    get_rhime_model,
    get_sampling_kwargs_from_rhime_outs,
    make_idata_from_rhime_outs,
)


def clean_rhime_output(ds: xr.Dataset) -> xr.Dataset:
    """Take raw RHIME output and rename/drop/create variables to get dataset ready for further processing."""
    ds = (
        ds.rename_vars(stepnum="draw", paramnum="nlatent", numBC="nBC", measurenum="nmeasure")
        .drop_dims(["nUI", "nlatent"])
        .swap_dims(nsite="nsites", steps="draw")
        .rename(
            {
                "nsites": "nsite",
                "nparam": "nx",
                "nBC": "nbc",
                "xtrace": "x",
                "bctrace": "bc",
                "sigtrace": "sigma",
            }
        )
    )
    ds["x"] = ds.x.assign_coords(nx=("nx", np.arange(ds.dims["nx"])))
    ds["mu_bc"] = (ds.bcsensitivity @ ds.bc).transpose("draw", ...)

    ds = ds[
        [
            "Yobs",
            "Yerror",
            "Ytime",
            "x",
            "bc",
            "sigma",
            "mu_bc",
            "siteindicator",
            "sigmafreqindex",
            "sitenames",
            "fluxapriori",
            "basisfunctions",
            "xsensitivity",
            "bcsensitivity",
        ]
    ]

    return ds


def nmeasure_to_site_time_data_array(
    da: xr.DataArray, site_indicators: xr.DataArray, site_names: xr.DataArray, times: xr.DataArray
) -> xr.DataArray:
    """Convert `nmeasure` dimension to multi-index over `site` and `time.`"""
    site_dict = dict(site_names.to_series())
    da = (
        xr.concat(
            [
                da.where(site_indicators == site_num, drop=True)
                .assign_coords(nmeasure=times.where(site_indicators == site_num, drop=True))
                .rename(nmeasure="time")
                .expand_dims({"site": [site_code]})
                for site_num, site_code in site_dict.items()
            ],
            dim="site",
        )
        .stack(nmeasure=["site", "time"])
        .dropna("nmeasure")
        .transpose("nmeasure", ...)
    )

    return da


def nmeasure_to_site_time(
    ds: xr.Dataset, site_indicators: xr.DataArray, site_names: xr.DataArray, times: xr.DataArray
) -> xr.Dataset:
    """Convert `nmeasure` dimension to multi-index over `site` and `time.`"""
    time_vars = [dv for dv in ds.data_vars if "nmeasure" in list(ds[dv].coords)]
    ds_tv = ds[time_vars]
    ds_no_tv = ds.drop_dims("nmeasure")

    site_dict = dict(site_names.to_series())
    ds_tv = (
        xr.concat(
            [
                ds_tv.where(site_indicators == site_num, drop=True)
                .expand_dims({"site": [site_code]})
                .assign_coords(nmeasure=times.where(site_indicators == site_num, drop=True))
                .rename_vars(nmeasure="time")
                for site_num, site_code in site_dict.items()
            ],
            dim="site",
        )
        .swap_dims(nmeasure="time")
        .stack(nmeasure=["site", "time"])
        .dropna("nmeasure")
        .transpose("nmeasure", ...)
    )

    return xr.merge([ds_no_tv, ds_tv])


InvOut = TypeVar("InvOut", bound="InversionOutput")  # for classmethod


@dataclass
class InversionOutput:
    """dataclass to hold the quantities we need to calculate outputs."""
    obs: xr.DataArray
    site_coordinates: xr.Dataset
    flux: xr.DataArray
    basis: xr.DataArray
    model: pm.Model
    trace: az.InferenceData
    site_indicators: xr.DataArray
    site_names: xr.DataArray
    times: xr.DataArray

    @classmethod
    def from_rhime(cls: type[InvOut], ds: xr.Dataset, min_model_error: float, ndraw: Optional[int] = None) -> InvOut:
        """Make InversionOutput object from RHIME output dataset."""
        flux = ds.fluxapriori
        site_coordinates = ds[["sitelons", "sitelats"]]

        ds_clean = clean_rhime_output(ds)
        site_indicators = ds_clean.siteindicator

        basis = get_xr_dummies(ds_clean.basisfunctions, cat_dim="nx")

        model_kwargs = get_sampling_kwargs_from_rhime_outs(ds)
        model = get_rhime_model(
            ds_clean,
            xprior_dims="nx",
            bcprior_dims="nbc",
            sigprior_dims=("nsigma_site", "nsigma_time"),
            min_model_error=min_model_error,
            **model_kwargs  # type: ignore
        )

        if ndraw is not None:
            trace = make_idata_from_rhime_outs(ds_clean, ndraw=ndraw)
        else:
            trace = make_idata_from_rhime_outs(ds_clean)

        return cls(
            obs=ds_clean.Yobs,
            site_coordinates=site_coordinates,
            flux=flux,
            basis=basis,
            model=model,
            trace=trace,
            site_indicators=site_indicators,
            site_names=ds_clean.sitenames,
            times=ds_clean.Ytime,
        )

    def sample_predictive_distributions(self, ndraw: int = 1000) -> None:
        """Sample prior and posterior predictive distributions.

        This creates prior samples as a side-effect.
        """
        self.trace.extend(pm.sample_prior_predictive(ndraw, self.model))
        self.trace.extend(pm.sample_posterior_predictive(self.trace, model=self.model, var_names=["y"]))

    def get_trace_dataset(
        self, convert_nmeasure: bool = True, var_names: Optional[Union[str, list[str]]] = None
    ) -> xr.Dataset:
        """Return an xarray Dataset containing a prior/posterior parameter/predictive samples.

        Args:
            convert_nmeasure: if True, convert `nmeasure` coordinate to multi-index comprising `time` and `site`.
            var_names: (list of) variables to select. For instance, "x" will return "x_prior" and "x_posterior".

        Returns:
            xarray Dataset containing a prior/posterior parameter/predictive samples.
        """
        trace_ds = convert_idata_to_dataset(self.trace)

        if convert_nmeasure:
            trace_ds = nmeasure_to_site_time(trace_ds, self.site_indicators, self.site_names, self.times)

        if var_names is not None:
            if isinstance(var_names, str):
                var_names = [var_names]

            data_vars = []
            for dv in trace_ds.data_vars:
                for name in var_names:
                    if str(dv).startswith(name):
                        data_vars.append(dv)

            trace_ds = trace_ds[data_vars]

        return trace_ds

    def start_time(self) -> np.datetime64:
        """Return start date of inversion."""
        return self.times.min().values  # type: ignore

    def period_midpoint(self) -> np.datetime64:
        """Return midpoint of inversion period."""
        half_of_period = (self.times.max().values - self.times.min().values) / 2
        return self.times.min().values + half_of_period  # type: ignore

    def get_obs(self, unstack_nmeasure: bool = True) -> xr.DataArray:
        """Return y observations.

        By default, `nmeasure` is converted to `site` and `time`.
        """
        result = nmeasure_to_site_time_data_array(self.obs, self.site_indicators, self.site_names, self.times)

        if unstack_nmeasure:
            return result.unstack("nmeasure")

        return result
