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


def get_rhime_model2(
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


def get_prior_samples(model: pm.Model) -> az.InferenceData:
    """Sample prior and prior predictive distributions from RHIME model."""
    idata = pm.sample_prior_predictive(1000, model=model)

    return cast(az.InferenceData, idata)


def make_idata_from_rhime_outs(rhime_out_ds: xr.Dataset, ndraw: int = 1000) -> az.InferenceData:
    """Create arviz InferenceData with posterior group created from RHIME output."""
    trace_dvs = [dv for dv in rhime_out_ds.data_vars if "draw" in list(rhime_out_ds[dv].coords)]
    traces = rhime_out_ds[trace_dvs].expand_dims({"chain": [0]})

    if traces.sizes["draw"] > ndraw:
        thin_by = traces.sizes["draw"] // ndraw
        traces = traces.isel(draw=slice(None, None, thin_by))
        traces = traces.assign_coords(draw=(traces.draw / thin_by).astype(int))

    return az.InferenceData(posterior=traces)


def get_posterior_predictive_samples(
    idata: az.InferenceData,
    model: pm.Model,
    var_names: list[str] = ["y", "muBC"],
    thin_by: Optional[int] = None,
) -> az.InferenceData:
    """Sample prior and prior predictive distributions from RHIME model.

    The default variables that are sampled are:
        "y", the modelled observations including modelling, and
        "muBC", the result of the forward model applied to the posterior
            boundary conditions trace (this does not include any model error).

    This will return as many samples as are in the InferenceData, which may take a fairly
    long time. The number of samples can be reduced by passing a value for `thin_by`.

    Args:
        idata: posterior trace from `model`.
        model: PyMC model (must have an observed variable).
        var_names: list of variables from `model` to sample from.
        thin_by: factor to thin the idata by; this will reduce the number of samples
            generated by the same factor.

    Returns:
        InferenceData with posterior predictive samples.
    """
    if thin_by is None:
        post_preds = pm.sample_posterior_predictive(idata, model=model, var_names=var_names)
    else:
        post_preds = pm.sample_posterior_predictive(
            idata.isel(draw=slice(None, None, thin_by)), model=model, var_names=var_names
        )

    return cast(az.InferenceData, post_preds)


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


def get_all_traces(
    outs: list[xr.Dataset], min_model_error: float, thin_by: Optional[int] = None
) -> list[xr.Dataset]:
    """Given a list of RHIME output datasets, construct corresponding datasets containing traces for prior,
    posterior, prior predictive, and posterior predictive distributions.
    """
    # reconstruct RHIME model for each RHIME output dataset
    model_kwargs = [get_sampling_kwargs_from_rhime_outs(ds, min_model_error=min_model_error) for ds in outs]
    models = [get_rhime_model(ds, **kwargs) for ds, kwargs in zip(outs, model_kwargs)]  # type: ignore

    # reconstruct posterior InferenceData from RHIME output datasets
    idatas = [make_idata_from_rhime_outs(ds) for ds in outs]

    # get prior and prior predictive samples
    for idata, model in zip(idatas, models):
        idata.extend(get_prior_samples(model))

    # get posterior predictive samples
    for idata, model in zip(idatas, models):
        idata.extend(get_posterior_predictive_samples(idata=idata, model=model, thin_by=thin_by))

    result = [convert_idata_to_dataset(idata) for idata in idatas]

    return result
