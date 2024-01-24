#!/usr/bin/env python
from typing import cast, Optional, Union, TypeVar

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
    function_dict = {cd.lower(): cd for cd in continuous.__all__}
    params = prior_params.copy()
    pdf = str(params.pop("pdf")).lower()  # str is just for typing...
    try:
        dist = getattr(continuous, function_dict[pdf])
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

    Args:
        rhime_outs_ds: output of RHIME inversion
        xprior: dictionary containing description of x prior pdf (flux uncertainty)
        bcprior: dictionary containing description of bc prior pdf (boundary conditions uncertainty)
        sigprior: dictionary containing description of sig prior pdf (observation uncertainty)
        min_model_error: constant minimum model error
        coord_dim: name of dimension used for coordinates of observations
        rename_coords: rename model coordinates so that samples created by this model have the same
            coordinates as `rhime_outs_ds`.
        x_name: name of xprior variable; this name is used for the corresponding data variable in the
            sampling output of this model.
        bc_name: name of bcprior variable; this name is used for the corresponding data variable in the
            sampling output of this model.
        sig_name: name of sigprior variable; this name is used for the corresponding data variable in the
            sampling output of this model.

    Returns:
        PyMC model for RHIME model with given priors and minimum model error.
    """
    if xprior is None:
        xprior = {"pdf": "truncatednormal", "mu": 1.0, "sigma": 1.0, "lower": 0.0}
    if bcprior is None:
        bcprior = {"pdf": "truncatednormal", "mu": 1.0, "sigma": 0.1, "lower": 0.0}
    if sigprior is None:
        sigprior = {"pdf": "uniform", "lower": 0.1, "upper": 1.0}

    coords_dict = {coord_dim: rhime_outs_ds[coord_dim]}
    if rename_coords:
        for dim in ["nparam", "nBC", "nsigma_site", "nsigma_time"]:
            coords_dict[dim] = rhime_outs_ds[dim]
        xprior_kwargs = {"dims": "nparam"}
        bcprior_kwargs = {"dims": "nBC"}
        sigprior_kwargs = {"dims": ("nsigma_site", "nsigma_time")}
    else:
        xprior_kwargs = {"shape": rhime_outs_ds.sizes["nparam"]}
        bcprior_kwargs = {"shape": rhime_outs_ds.sizes["nBC"]}
        sigprior_kwargs = {"shape": (rhime_outs_ds.sizes["nsigma_site"], rhime_outs_ds.sizes["nsigma_time"])}

    with pm.Model(coords=coords_dict) as model:
        x = parse_prior(x_name, xprior, **xprior_kwargs)
        bc = parse_prior(bc_name, bcprior, **bcprior_kwargs)
        sigma = parse_prior(sigma_name, sigprior, **sigprior_kwargs)

        muBC = pm.Deterministic("muBC", pt.dot(rhime_outs_ds.bcsensitivity.values, bc), dims=coord_dim)
        mu = pt.dot(rhime_outs_ds.xsensitivity.values, x) + muBC

        mult_error = (
            np.abs(pt.dot(rhime_outs_ds.xsensitivity.values, x))
            * sigma[
                rhime_outs_ds.siteindicator.values.astype(int),
                rhime_outs_ds.sigmafreqindex.values.astype(int),
            ]
        )
        epsilon = pt.sqrt(rhime_outs_ds.Yerror.values**2 + mult_error**2 + min_model_error**2)
        y = pm.Normal("y", mu, epsilon, dims=coord_dim, observed=rhime_outs_ds.Yobs)

    return model


def get_prior_samples(
    rhime_out_ds: Optional[xr.Dataset] = None,
    model: Optional[pm.Model] = None,
    xprior: Optional[PriorArgs] = None,
    bcprior: Optional[PriorArgs] = None,
    sigprior: Optional[PriorArgs] = None,
    min_model_error: float = 20.0,
    coord_dim: str = "nmeasure",
) -> az.InferenceData:
    """Sample prior and prior predictive distributions from RHIME model."""
    if model is None:
        if rhime_out_ds is None:
            raise ValueError("If `model` is not provided, then `rhime_out_ds` must be provided.")
        model = get_rhime_model(rhime_out_ds, xprior, bcprior, sigprior, min_model_error, coord_dim)

    idata = pm.sample_prior_predictive(1000, model=model)

    return cast(az.InferenceData, idata)


def make_idata_from_rhime_outs(rhime_out_ds: xr.Dataset) -> az.InferenceData:
    """Create arviz InferenceData with posterior group created from RHIME output."""
    trace_dvs = [dv for dv in rhime_out_ds.data_vars if "trace" in str(dv)]
    return az.InferenceData(posterior=rhime_out_ds[trace_dvs].expand_dims({"chain": [0]}))


def get_posterior_predictive_samples(
    rhime_out_ds: Optional[xr.Dataset] = None,
    model: Optional[pm.Model] = None,
    idata: Optional[az.InferenceData] = None,
    xprior: Optional[PriorArgs] = None,
    bcprior: Optional[PriorArgs] = None,
    sigprior: Optional[PriorArgs] = None,
    min_model_error: float = 20.0,
    coord_dim: str = "nmeasure",
) -> az.InferenceData:
    """Sample prior and prior predictive distributions from RHIME model."""
    if model is None or idata is None:
        if rhime_out_ds is None:
            raise ValueError(
                "If `rhime_out_ds` is not provided, then both `model` and `idata` must be provided."
            )
        else:
            if model is None:
                model = get_rhime_model(rhime_out_ds, xprior, bcprior, sigprior, min_model_error, coord_dim)

            if idata is None:
                idata = make_idata_from_rhime_outs(rhime_out_ds)

    post_preds = cast(
        az.InferenceData, pm.sample_posterior_predictive(idata, model=model, var_names=["y", "muBC"])
    )

    return post_preds


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
