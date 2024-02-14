#!/usr/bin/env python
from typing import Optional, Union, TypeVar

import arviz as az
import numpy as np
import pymc as pm
import pytensor.tensor as pt
import xarray as xr
from pymc.distributions import continuous
from pytensor.tensor.variable import TensorVariable
from xarray.core.common import DataWithCoords


# type for xr.Dataset *or* xr.DataArray
xrData = TypeVar("xrData", bound=DataWithCoords)

# type alias for prior args
PriorArgs = dict[str, Union[str, float]]


def parse_prior(name: str, prior_params: PriorArgs, **kwargs) -> TensorVariable:
    """
    Parses all PyMC continuous distributions:
    https://docs.pymc.io/api/distributions/continuous.html

    Args:
        name:
          name of variable in the pymc model
        prior_params:
          dict of parameters for the distribution, including 'pdf' for the distribution to use.
          The value of `prior_params["pdf"]` must match the name of a PyMC continuous
          distribution: https://docs.pymc.io/api/distributions/continuous.html
        **kwargs: for instance, `shape` or `dims`

    Returns:
        continuous PyMC distribution

    For example:
    ```
    params = {"pdf": "uniform", "lower": 0.0, "upper": 1.0}
    parse_prior("x", params, shape=(20, 20))
    ```
    will create a 20 x 20 array of uniform random variables.

    Alternatively,
    ```
    params = {"pdf": "uniform", "lower": 0.0, "upper": 1.0}
    parse_prior("x", params, dims="nmeasure"))
    ```
    will create an array of uniform random variables with the same shape
    as the dimension coordinate `nmeasure`. This can be used if `pm.Model`
    is provided with coordinates.

    Note: calling `parse_prior` inside a `pm.Model` context (i.e. after `with pm.Model()`)
    has an important side-effect of registering the random variable with the model. Typically,
    `parse_prior` will be used in a `pm.Model` context, although it could be used to construct
    random variables for other purposes.
    """
    # create dict to lookup continuous PyMC distributions by name, ignoring case
    pdf_dict = {cd.lower(): cd for cd in continuous.__all__}

    params = prior_params.copy()
    pdf = str(params.pop("pdf")).lower()  # str is just for typing...
    try:
        dist = getattr(continuous, pdf_dict[pdf])
    except AttributeError:
        raise ValueError(
            f"The distribution '{pdf}' doesn't appear to be a continuous distribution defined by PyMC."
        )
    return dist(name, **params, **kwargs)


def get_sampling_kwargs_from_rhime_outs(
    rhime_outs: xr.Dataset, min_model_error: Optional[float] = None
) -> dict[str, Union[float, PriorArgs]]:
    """Make dict of arguments to use with `get_rhime_model` given a RHIME output file
    and minimum model error.

    Convenience function for creating PARIS outputs based on RHIME/HBMCMC outputs.
    """

    def get_prior_from_attrs(attr_val: str) -> PriorArgs:
        """Get dict of prior arguments from RHIME output attributes."""
        prior = {}
        split = attr_val.split(",")
        for k, v in zip(split[::2], split[1::2]):
            prior[k] = v if k == "pdf" else float(v)
        return prior

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
    rhime_outs_ds: xr.Dataset,
    xprior: PriorArgs,
    bcprior: PriorArgs,
    sigprior: PriorArgs,
    xprior_dims: Union[str, tuple[str, ...]],
    bcprior_dims: Union[str, tuple[str, ...]],
    sigprior_dims: Union[str, tuple[str, ...]],
    min_model_error: float = 20.0,
    coord_dim: str = "nmeasure",
) -> pm.Model:
    """Make RHIME model wih model error given by multiplicative scaling of pollution events
    plus constant minimum value.

    Args:
        rhime_outs_ds: output of RHIME inversion
        xprior: dictionary containing description of x prior pdf (flux uncertainty)
        bcprior: dictionary containing description of bc prior pdf (boundary conditions uncertainty)
        sigprior: dictionary containing description of sig prior pdf (observation uncertainty)
        xprior_dims: coordinate dims from rhime_outs_ds to use
        bcprior_dims: coordinate dims from rhime_outs_ds to use
        sigprior_dims: coordinate dims from rhime_outs_ds to use
        min_model_error: constant minimum model error
        coord_dim: name of dimension used for coordinates of observations

    Returns:
        PyMC model for RHIME model with given priors and minimum model error.
    """
    coords_dict = {coord_dim: rhime_outs_ds[coord_dim]}
    for dim in [xprior_dims, bcprior_dims, sigprior_dims]:
        if isinstance(dim, str):
            coords_dict[dim] = rhime_outs_ds[dim]
        else:
            for d in dim:
                coords_dict[d] = rhime_outs_ds[d]

    with pm.Model(coords=coords_dict) as model:
        x = parse_prior("x", xprior, dims=xprior_dims)
        bc = parse_prior("bc", bcprior, dims=bcprior_dims)
        sigma = parse_prior("sigma", sigprior, dims=sigprior_dims)

        mu_bc = pm.Deterministic("mu_bc", pt.dot(rhime_outs_ds.bcsensitivity.values, bc), dims=coord_dim)
        mu = pt.dot(rhime_outs_ds.xsensitivity.values, x) + mu_bc

        mult_error = (
            np.abs(pt.dot(rhime_outs_ds.xsensitivity.values, x))
            * sigma[
                rhime_outs_ds.siteindicator.values.astype(int),
                rhime_outs_ds.sigmafreqindex.values.astype(int),
            ]
        )
        epsilon = pm.Deterministic(
            "epsilon",
            pt.sqrt(rhime_outs_ds.Yerror.values**2 + mult_error**2 + min_model_error**2),
            dims=coord_dim,
        )
        pm.Normal("y", mu, epsilon, dims=coord_dim, observed=rhime_outs_ds.Yobs)

    return model


def make_idata_from_rhime_outs(rhime_out_ds: xr.Dataset, ndraw: int = 1000) -> az.InferenceData:
    """Create arviz InferenceData with posterior group created from RHIME output."""
    trace_dvs = [dv for dv in rhime_out_ds.data_vars if "draw" in list(rhime_out_ds[dv].coords)]
    traces = rhime_out_ds[trace_dvs].expand_dims({"chain": [0]})

    if traces.sizes["draw"] > ndraw:
        thin_by = traces.sizes["draw"] // ndraw
        traces = traces.isel(draw=slice(None, None, thin_by))
        traces = traces.assign_coords(draw=(traces.draw / thin_by).astype(int))

    return az.InferenceData(posterior=traces)


def convert_idata_to_dataset(idata: az.InferenceData) -> xr.Dataset:
    """Merge prior, prior predictive, posterior, and posterior predictive samples into a single
    xr.Dataset.
    """
    traces = []
    for group in idata.groups():
        if "prior" in group or "posterior" in group:
            trace = idata[group]
            rename_dict = {dv: f"{dv}_{group}" for dv in trace.data_vars}
            traces.append(trace.rename_vars(rename_dict).squeeze("chain", drop=True))
    return xr.merge(traces)
