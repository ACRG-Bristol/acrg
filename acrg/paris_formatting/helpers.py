#!/usr/bin/env python
import json
from typing import Any, cast, Optional, Sequence, Union, TypeVar

import arviz as az
import numpy as np
import openghg_inversions as oi
import pandas as pd
import pymc as pm
import sparse
import xarray as xr
from openghg_inversions import convert, utils
from scipy import stats
from sparse import COO
from xarray.core.common import DataWithCoords


# type for xr.Dataset *or* xr.DataArray
xrData = TypeVar("xrData", bound=DataWithCoords)


def make_quantiles(
    da: xr.DataArray, probs: Sequence[float] = [0.025, 0.159, 0.841, 0.975], sample_dim="steps"
) -> xr.DataArray:
    """Return xr.DataArray of quantiles computed over dimension `sample_dim`.

    NOTE: xarray has a .quantiles method built in, maybe we should use this.
    """
    probs_da = xr.DataArray(probs, coords=[probs], dims=["probs"])

    # make function to apply
    # we will pass `q=probs_da` so that the coordinates of probs will be propegated
    def func(a, q):
        qs = np.quantile(a, q, axis=-1)  # apply along input_core_dim = sample_dim
        qs = qs[..., 0]  # contracted dimension at axis=-1 is left with length 1, need to remove it
        qs = np.moveaxis(qs, 0, -1)
        return qs

    result = xr.apply_ufunc(func, da, probs_da, input_core_dims=[[sample_dim], []])
    return result.transpose("probs", ...)  # we want "probs" first


def calc_mode(da: xr.DataArray, sample_dim="steps") -> xr.DataArray:
    """Calculate the (KDE smoothed) mode of a data array containing MCMC
    samples.

    This can be parallelized if you chunk the DataArray first, e.g.

    da_chunked = da.chunk({"basis_region": 10})
    """

    def mode_of_row(row):
        if np.all(np.isnan(row)):
            return np.nan
        elif np.nanmax(row) > np.nanmin(row):
            xvals = np.linspace(np.nanmin(row), np.nanmax(row), 200)
            kde = stats.gaussian_kde(row).evaluate(xvals)
            return xvals[kde.argmax()]
        else:
            return np.nanmean(row)

    def func(arr):
        return np.apply_along_axis(func1d=mode_of_row, axis=-1, arr=arr)

    return xr.apply_ufunc(func, da, input_core_dims=[[sample_dim]], dask="parallelized")


def get_xr_dummies(
    da: xr.DataArray,
    categories: Optional[Union[Sequence[Any], pd.Index, xr.DataArray, np.ndarray]] = None,
    cat_dim: str = "categories",
    sparse: bool = True,
):
    """Create 0-1 dummy matrix from DataArray with values that correspond to categories.

    If the values of `da` are integers 0-N, then the result has N + 1 columns, and the (i, j) coordiante
    of the result is 1 if `da[i] == j`, and is 0 otherwise.

    Args:
        da: DataArray encoding categories.
        categories: optional coordinates for categories.
        cat_dim: dimension for categories coordinate
        sparse: if True, store values in sparse.COO matrix

    Returns:
        Dummy matrix corresponding to the input vector. Its dimensions are the same as the
    input DataArray, plus an additional "categories" dimension, which  has one value for each
    distinct value in the input DataArray.
    """
    # stack if `da` is not one dimensional
    stack_dim = ""
    if len(da.dims) > 1:
        stack_dim = "".join([str(dim) for dim in da.dims])
        da = da.stack({stack_dim: da.dims})

    dummies = pd.get_dummies(da.values, dtype=int, sparse=sparse)

    # put dummies into DataArray with the right coords and dims
    values = COO.from_scipy_sparse(dummies.sparse.to_coo()) if sparse else dummies.values
    if categories is None:
        categories = np.arange(values.shape[1])
    coords = da.coords.merge({cat_dim: categories}).coords  # coords.merge returns Dataset, we want the coords
    result = xr.DataArray(values, coords=coords)

    # if we stacked `da`, unstack result before returning
    return result.unstack(stack_dim) if stack_dim else result


def sparse_xr_dot(
    da1: xr.DataArray, da2: xrData, debug: bool = False, broadcast_dims: Optional[Sequence[str]] = None
) -> xrData:
    """Compute the matrix "dot" of a tuple of DataArrays with sparse.COO values.

    This multiplies and sums over all common dimensions of the input DataArrays, and
    preserves the coordinates and dimensions that are not summed over.

    Common dimensions are automatically selected by name. The input arrays must  have at
    least one dimension in common. All matching dimensions will be used for multiplication.

    NOTE: this function shouldn't be necessary, but `da1 @ da2` doesn't work properly if the
    values of `da1` and `da2` are `sparse.COO` arrays.

    Args:
        da1, da2: xr.DataArrays to multiply and sum along common dimensions.
        debug: if true, will print the dimensions of the inputs to `sparse.tensordot`
            as well as the dimension of the result.
        along_dim: name

    Returns:
        xr.DataArray containing the result of matrix/tensor multiplication

    Raises:
        ValueError if the input DataArrays have no common dimensions to multiply.
    """
    common_dims = set(da1.dims).intersection(set(da2.dims))
    nc = len(common_dims)

    if nc == 0:
        raise ValueError(f"DataArrays \n{da1}\n{da2}\n have no common dimensions. Cannot compute `dot`.")

    if broadcast_dims is not None:
        _broadcast_dims = set(broadcast_dims).intersection(common_dims)
    else:
        _broadcast_dims = set([])

    contract_dims = common_dims.difference(_broadcast_dims)
    ncontract = len(contract_dims)

    tensor_dot_axes = tuple([tuple(range(-ncontract, 0))] * 2)
    input_core_dims = [list(contract_dims)] * 2

    # compute tensor dot on last nc coordinates (because core dims are moved to end)
    # and then drop 1D coordinates resulting from summing
    def _func(x, y, debug=False):
        result = sparse.tensordot(x, y, axes=tensor_dot_axes)  # type: ignore

        xs = list(x.shape[:-ncontract])
        nxs = len(xs)
        ys = list(y.shape[:-ncontract])
        nys = len(ys)
        pad_y = nxs - nys
        if pad_y > 0:
            ys = [1] * pad_y + ys

        idx1, idx2 = [], []
        for i, j in zip(xs, ys):
            if i == j or j == 1:
                idx1.append(slice(None))
                idx2.append(0)
            elif i == 1:
                # x broadcasted to match y's dim
                idx1.append(0)
                idx2.append(slice(None))

        if debug:
            print("pad y", pad_y)
            print(xs, ys)
            print(idx1, idx2)

        idx2 = idx2[pad_y:]
        idx3 = [0] * (result.ndim - len(idx1) - len(idx2))
        idx = tuple(idx1 + idx3 + idx2)

        if debug:
            print("x.shape", x.shape, "y.shape", y.shape)
            print("idx", idx)
            print("result shape:", result.shape)

        return result[idx]  # type: ignore

    def wrapper(da1, da2):
        for arr in [da1, da2]:
            print(f"_func received array of type {type(arr)}, shape {arr.shape}")
        result = _func(da1, da2, debug=True)
        print(f"_func result shape: {result.shape}\n")
        return result

    func = wrapper if debug else _func

    return xr.apply_ufunc(func, da1, da2, input_core_dims=input_core_dims, join="outer")


def get_area_grid(lat: xr.DataArray, lon: xr.DataArray) -> xr.DataArray:
    """Return xr.DataArray with coordinate dimensions ("lat", "lon") containing
    the area of each grid cell centered on the coordinates.

    Args:
        lat: latitude values
        lon: longitude values

    Returns:
        xr.DataArray with grid cell areas.
    """
    ag_vals = utils.areagrid(lat.values, lon.values)
    return xr.DataArray(ag_vals, coords=[lat, lon], dims=["lat", "lon"], name="area_grid")


def get_x_to_country_mat(
    countries: xr.Dataset,
    hbmcmc_outs: Optional[xr.Dataset] = None,
    flux: Optional[xr.DataArray] = None,
    area_grid: Optional[xr.DataArray] = None,
    basis_functions: Optional[xr.DataArray] = None,
    sparse: bool = False,
    basis_cat_dim: str = "basis_region",
    country_selection: Optional[list[str]] = None,
) -> xr.DataArray:
    """Construct a sparse matrix mapping from x sensitivities to country totals.

    Args:
        countries: xr.Dataset from country file. Must have variables: "country" with coordinate
            dimensions ("lat", "lon"), and "name" (with no coordinate dimensions).
        hbmcmc_outs: xr.Dataset from `hbmcmc_postprocessouts`.
        flux: flux used in inversion. If a constant flux was used, then `aprioriflux` from
            `hbmcmc_outs` can be used.
        area_grid: areas of each grid cell in inversion domain.
        basis_functions: xr.DataArray with coordinate dimensions ("lat", "lon") whose values assign
            grid cells to basis function boxes.
        sparse: if True, values of returned DataArray are `sparse.COO` array.

    Returns:
        xr.DataArray with coordinate dimensions ("country", "basis_region")
    """
    # we need flux, area_grid, and basis_functions. Get them from hbmcmc_outs if not provided.
    if any(arg is None for arg in [flux, area_grid, basis_functions]):
        if hbmcmc_outs is None:
            raise ValueError(
                "If `hbmcmc_outs` is given, you must provide `flux`, `area_grid`, and `basis_functions`"
            )
        else:
            if flux is None:
                flux = hbmcmc_outs.fluxapriori
            if area_grid is None:
                area_grid = get_area_grid(hbmcmc_outs.lat, hbmcmc_outs.lon)
            if basis_functions is None:
                basis_functions = hbmcmc_outs.basisfunctions

    # at this point, flux, area_grid, and basis_functions are not None
    flux = cast(xr.DataArray, flux)
    area_grid = cast(xr.DataArray, area_grid)
    basis_functions = cast(xr.DataArray, basis_functions)

    # create dummy matrices from country  DataArrays
    country_mat = get_xr_dummies(countries.country, cat_dim="country", categories=countries.name)

    # filter based on countries (better performance by doing this early)
    if country_selection:
        filt = countries.name.isin(country_selection)
        country_mat = country_mat.where(filt, drop=True)

    # create dummy matrices from basis DataArrays
    basis_mat = get_xr_dummies(basis_functions, cat_dim=basis_cat_dim)

    # compute matrix/tensor product: country_mat.T @ (area_grid * flux * basis_mat)
    # transpose doesn't need to be taken explicitly because alignment is done by dimension name
    result = sparse_xr_dot(country_mat, area_grid * flux * basis_mat)

    if sparse:
        return result

    # hack since `.to_numpy()` doesn't work right with sparse arrays
    return xr.apply_ufunc(lambda x: x.todense(), result)


def get_country_trace(
    species: str,
    x_trace: xrData,
    x_to_country: xr.DataArray,
) -> xrData:
    """Return trace for total country emissions.

    Args:
        species: name of species, e.g. "co2", "ch4", "sf6", etc.
        x_trace: xr.DataArray or xr.Dataset with coordinate dimensions ("draw", "nparam")
        x_to_country: xr.DataArray with result from `get_x_to_country_mat`

    Returns:
        xr.DataArray with coordinates ("country", "stepnum") and dimensions ("country", "steps")

    TODO: add attributes?
    TODO: there is a "country unit" conversion in the old code, but it seems to always product
          1.0, based on how it is used in hbmcmc
    TODO: PyMC uses "draw" (or "draws"?) instead of "steps", so we might need to handle both or
          change "hbmcmc_postprocessouts". Ideally, this trace would be compatible with arviz.InferenceData

    DONE: make a version where `x_to_country` can be passed as an argument...
    NOTE: if number of basis regions changes then this is no good...
    """
    raw_trace = sparse_xr_dot(x_to_country, x_trace)
    molar_mass = convert.molar_mass(species)
    return raw_trace * 365 * 24 * 3600 * molar_mass


def convert_time_to_unix_epoch(x: xrData) -> xrData:
    """Convert `time` coordinate of xarray Dataset or DataArray to number of seconds since
    1 Jan 1970 (the "UNIX epoch").
    """
    return x.assign_coords(time=(pd.DatetimeIndex(x.time) - pd.Timestamp("1970-01-01")) // pd.Timedelta("1s"))


def convert_unix_epoch_to_time(x: xrData) -> xrData:
    """Convert `time` coordinate of xarray Dataset or DataArray to number of seconds since
    1 Jan 1970 (the "UNIX epoch").
    """
    return x.assign_coords(time=(pd.Timedelta("1s") * x.time + pd.Timestamp("1970-01-01")))
