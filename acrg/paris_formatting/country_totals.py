#!/usr/bin/env python3
import json
from typing import Any, cast, Optional, Sequence, Union

import numpy as np
import openghg_inversions as oi
from openghg_inversions import utils
import pandas as pd
import sparse
from sparse import COO
import xarray as xr


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


def sparse_xr_dot(da1: xr.DataArray, da2: xr.DataArray, debug: bool = False) -> xr.DataArray:
    """Compute the matrix "dot" of a tuple of DataArrays with sparse.COO values.

    This multiplies and sums over all common dimensions of the input DataArrays, and
    preserves the coordinates and dimensions that are not summed over.

    Common dimensions are automatically selected by name. The input arrays must  have at
    least one dimension in common. All matching dimensions will be used for multiplication.

    NOTE: this function shouldn't be necessary, but `da1 @ da2` doesn't work properly if the
    values of `da1` and `da2` are `sparse.COO` arrays.

    Args:
        *arrays: xr.DataArrays to multiply and sum along common dimensions.
        debug: if true, will print the dimensions of the inputs to `sparse.tensordot`
            as well as the dimension of the result.

    Returns:
        xr.DataArray containing the result of matrix/tensor multiplication

    Raises:
        ValueError if the input DataArrays have no common dimensions to multiply.
    """
    common_dims = set(da1.dims).intersection(set(da2.dims))
    nc = len(common_dims)

    if nc == 0:
        raise ValueError(f"DataArrays \n{da1}\n{da2}\n have no common dimensions. Cannot compute `dot`.")

    tensor_dot_axes = tuple([tuple(range(-nc, 0))] * 2)
    input_core_dims = [list(common_dims)] * 2

    # compute tensor dot on last nc coordinates (because core dims are moved to end)
    # and then drop 1D coordinates resulting from summing
    def _func(x, y):
        result = sparse.tensordot(x, y, axes=tensor_dot_axes)  # type: ignore

        n1, n2, nr = len(da1.dims) - nc, len(da2.dims) - nc, len(result.shape)
        # we should have n1 axes from da1, n2 axes from da2, and nr - n1 - n2
        # axes left over (sometimes axes of length 1 are left by tensordot)
        idx = tuple([slice(None)] * n1 + [0] * (nr - n1 - n2) + [slice(None)] * n2)
        return result[idx]  # type: ignore

    def wrapper(da1, da2):
        for arr in [da1, da2]:
            print(f"_func received array of type {type(arr)}, shape {arr.shape}")
        result = _func(da1, da2)
        print(f"_func result shape: {result.shape}")
        return result

    func = wrapper if debug else _func

    return xr.apply_ufunc(func, da1, da2, input_core_dims=input_core_dims)


def get_area_grid_data_array(lat: xr.DataArray, lon: xr.DataArray) -> xr.DataArray:
    """Return xr.DataArray with coordinate dimensions ("lat", "lon") containing
    the area of each grid cell centered on the coordinates.

    Args:
        lat: latitude values
        lon: longitude values

    Returns:
        xr.DataArray with grid cell areas.
    """
    ag_vals = oi.utils.areagrid(lat.values, lon.values)
    return xr.DataArray(ag_vals, coords=[lat, lon], dims=["lat", "lon"], name="area_grid")


def get_x_to_country_mat(
    countries: xr.Dataset,
    hbmcmc_outs: Optional[xr.Dataset] = None,
    flux: Optional[xr.DataArray] = None,
    area_grid: Optional[xr.DataArray] = None,
    basis_functions: Optional[xr.DataArray] = None,
    sparse: bool = False,
):
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
                area_grid = get_area_grid_data_array(hbmcmc_outs.lat, hbmcmc_outs.lon)
            if basis_functions is None:
                basis_functions = hbmcmc_outs.basis_functions

    # at this point, flux, area_grid, and basis_functions are not None
    flux = cast(xr.DataArray, flux)
    area_grid = cast(xr.DataArray, area_grid)
    basis_functions = cast(xr.DataArray, basis_functions)

    # create dummy matrices from country and basis DataArrays
    country_mat = get_xr_dummies(countries.country, cat_dim="country", categories=countries.name)
    basis_mat = get_xr_dummies(basis_functions, cat_dim="basis_region")

    # compute matrix/tensor product: country_mat.T @ (area_grid * flux * basis_mat)
    # transpose doesn't need to be taken explicitly because alignment is done by dimension name
    result = sparse_xr_dot(country_mat, area_grid * flux * basis_mat)

    if sparse:
        return result

    # hack since `.to_numpy()` doesn't work right with sparse arrays
    return xr.apply_ufunc(lambda x: x.todense(), result)


def get_country_trace(
    countries: xr.Dataset,
    species: str,
    hbmcmc_outs: Optional[xr.Dataset] = None,
    flux: Optional[xr.DataArray] = None,
    area_grid: Optional[xr.DataArray] = None,
    basis_functions: Optional[xr.DataArray] = None,
    x_trace: Optional[xr.DataArray] = None,
):
    """Return trace for total country emissions.

    Args:
        countries: xr.Dataset from country file. Must have variables: "country" with coordinate
            dimensions ("lat", "lon"), and "name" (with no coordinate dimensions).
        hbmcmc_outs: xr.Dataset from `hbmcmc_postprocessouts`.
        flux: flux used in inversion. If a constant flux was used, then `aprioriflux` from
            `hbmcmc_outs` can be used.
        area_grid: areas of each grid cell in inversion domain.
        basis_functions: xr.DataArray with coordinate dimensions ("lat", "lon") whose values assign
            grid cells to basis function boxes.
        x_trace: xr.DataArray with coordinate dimensions ("steps", <other>), where <other> will
            be converted to "basis_region"

    Returns:
        xr.DataArray with coordinates ("country", "stepnum") and dimensions ("country", "steps")

    TODO: add attributes?
    TODO: there is a "country unit" conversion in the old code, but it seems to always product
          1.0, based on how it is used in hbmcmc
    TODO: PyMC uses "draw" (or "draws"?) instead of "steps", so we might need to handle both or
          change "hbmcmc_postprocessouts". Ideally, this trace would be compatible with arviz.InferenceData
    TODO: make a version where `x_to_country` can be passed as an argument...
    """
    if x_trace is None:
        if hbmcmc_outs is None:
            raise ValueError("If `hbmcmc_outs` is not given, `x_trace` must be provided.")
        else:
            x_trace = hbmcmc_outs.xtrace

    x_trace = cast(xr.DataArray, x_trace)

    # convert "other" coordinate of x_trace
    dim1, dim2 = x_trace.dims
    dim = dim2 if dim1 == "steps" or dim1 == "draw" else dim1
    x_trace = x_trace.rename_vars({dim: "basis_region"})

    x_to_country = get_x_to_country_mat(countries, hbmcmc_outs, flux, area_grid, basis_functions, sparse=True)

    raw_trace = sparse_xr_dot(x_to_country, x_trace)
    molar_mass = oi.convert.molar_mass(species)
    return raw_trace * 365 * 24 * 3600 * molar_mass
