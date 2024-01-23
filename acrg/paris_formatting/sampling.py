#!/usr/bin/env python
import json
from typing import Any, cast, Optional, Sequence, Union, TypeVar

import arviz as az
import numpy as np
import openghg_inversions as oi
import pandas as pd
import pymc as pm
import pytensor.tensor as pt
import sparse
import xarray as xr
from pymc.distributions import continuous
from openghg_inversions import convert
from openghg_inversions import utils
from scipy import stats
from sparse import COO
from xarray.core.common import DataWithCoords


# type for xr.Dataset *or* xr.DataArray
xrData = TypeVar("xrData", bound=DataWithCoords)

# type alias for prior args
PriorArgs = dict[str, Union[str, float]]


def parse_prior(name: str, prior_params: PriorArgs, **kwargs):
    """
    Parses all PyMC continuous distributions:
    https://docs.pymc.io/api/distributions/continuous.html

    Args:
        name:
          name of variable in the pymc model
        prior_params:
          dict of parameters for the distribution,
          including 'pdf' for the distribution to use
        **kwargs: for instance, `shape` or `dims`

    """
    function_dict = {cd.lower(): cd for cd in continuous.__all__}
    params = prior_params.copy()
    pdf = str(params.pop("pdf")).lower()  # str is just for typing...
    dist = getattr(continuous, function_dict[pdf])
    return dist(name, **params, **kwargs)


def get_prior_from_attrs(attr_val: str) -> PriorArgs:
    """Get dict of prior arguments from RHIME output attributes."""
    prior = {}
    split = attr_val.split(",")
    for k, v in zip(split[::2], split[1::2]):
        prior[k] = v if k == "pdf" else float(v)
    return prior


def get_sampling_kwargs(
    rhime_outs: xr.Dataset, min_model_error: Optional[float] = None
) -> dict[str, Union[float, PriorArgs]]:
    """Make dict of arguments to use with `get_rhime_model` given a RHIME output file
    and
    """
    result: dict[str, Union[float, PriorArgs]]
    result = dict(
        xprior=get_prior_from_attrs(rhime_outs.attrs["Emissions Prior"]),
        bcprior=get_prior_from_attrs(rhime_outs.attrs["BCs Prior"]),
        sigprior=get_prior_from_attrs(rhime_outs.attrs["Model error Prior"]),
    )
    if min_model_error:
        result["min_model_error"] = min_model_error
    return result


def get_rhime_model(
    ds: xr.Dataset,
    xprior: Optional[PriorArgs] = None,
    bcprior: Optional[PriorArgs] = None,
    sigprior: Optional[PriorArgs] = None,
    min_model_error: float = 20.0,
    coord_dim: str = "nmeasure",
    rename_coords: bool = True,
    x_name: str = "x",
    bc_name: str = "xbc",
    sigma_name: str = "sig",
) -> pm.Model:
    """Make RHIME model wih model error given by multiplicative scaling of pollution events
    plus constant minimum value.
    """
    if xprior is None:
        xprior = {"pdf": "truncatednormal", "mu": 1.0, "sigma": 1.0, "lower": 0.0}
    if bcprior is None:
        bcprior = {"pdf": "truncatednormal", "mu": 1.0, "sigma": 0.1, "lower": 0.0}
    if sigprior is None:
        sigprior = {"pdf": "uniform", "lower": 0.1, "upper": 1.0}

    coords_dict = {coord_dim: ds[coord_dim]}
    if rename_coords:
        for dim in ["nparam", "nBC", "nsigma_site", "nsigma_time"]:
            coords_dict[dim] = ds[dim]

    with pm.Model(coords=coords_dict) as model:
        x = parse_prior(x_name, xprior, dims="nparam")
        bc = parse_prior(bc_name, bcprior, dims="nBC")
        sigma = parse_prior(sigma_name, sigprior, dims=("nsigma_site", "nsigma_time"))

        muBC = pt.dot(ds.bcsensitivity.values, bc)
        mu = pt.dot(ds.xsensitivity.values, x) + muBC

        mult_error = (
            np.abs(pt.dot(ds.xsensitivity.values, x))
            * sigma[ds.siteindicator.values.astype(int), ds.sigmafreqindex.values.astype(int)]
        )
        epsilon = pt.sqrt(ds.Yerror.values**2 + mult_error**2 + min_model_error**2)
        ymodbc = pm.Normal("ymodbc", muBC, epsilon, dims=coord_dim)
        ymod = pm.Normal("ymod", mu, epsilon, dims=coord_dim)

    return model


def get_prior_samples(
    ds: xr.Dataset,
    xprior: Optional[PriorArgs] = None,
    bcprior: Optional[PriorArgs] = None,
    sigprior: Optional[PriorArgs] = None,
    min_model_error: float = 20.0,
    coord_dim: str = "nmeasure",
) -> az.InferenceData:
    """Sample prior and prior predictive distributions from RHIME model."""
    model = get_rhime_model(ds, xprior, bcprior, sigprior, min_model_error, coord_dim)
    idata = pm.sample_prior_predictive(1000, model=model)

    return cast(az.InferenceData, idata)


def combine_outs_and_prior_samples(
    outs: xr.Dataset, prior_samples: Union[xr.Dataset, az.InferenceData]
) -> xr.Dataset:
    """Combine trace variables from RHIME outputs and prior samples."""
    if isinstance(prior_samples, az.InferenceData):
        try:
            prior_samples = prior_samples.prior
        except AttributeError:
            raise ValueError("If `prior_samples` is arviz InferenceData, it must have a `prior` group.")
        else:
            prior_samples = cast(xr.Dataset, prior_samples).squeeze(drop=True)

    trace_dvs = [dv for dv in outs.data_vars if "trace" in str(dv)]
    prior_rename_dict = {dv: str(dv) + "_prior" for dv in prior_samples.data_vars}
    return xr.merge([outs[trace_dvs], prior_samples.rename_vars(prior_rename_dict)])
