"""
General methods for xarray Datasets and DataArrays.

`get_xr_dummies` applies pandas `get_dummies` to xarray DataArrays.

`sparse_xr_dot` multiplies a Dataset or DataArray by a DataArray
 with sparse underlying array. The built-in xarray functionality doesn't
work correctly.
"""
from typing import Any, Optional, Sequence, Union, TypeVar

import numpy as np
import pandas as pd
import sparse
import xarray as xr
from sparse import COO
from xarray.core.common import DataWithCoords


# type for xr.Dataset *or* xr.DataArray
DataSetOrArray = TypeVar("DataSetOrArray", bound=DataWithCoords)


def get_xr_dummies(
    da: xr.DataArray,
    categories: Optional[Union[Sequence[Any], pd.Index, xr.DataArray, np.ndarray]] = None,
    cat_dim: str = "categories",
    return_sparse: bool = True,
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

    dummies = pd.get_dummies(da.values, dtype=int, sparse=return_sparse)

    # put dummies into DataArray with the right coords and dims
    values = COO.from_scipy_sparse(dummies.sparse.to_coo()) if return_sparse else dummies.values
    if categories is None:
        categories = np.arange(values.shape[1])
    coords = da.coords.merge({cat_dim: categories}).coords  # coords.merge returns Dataset, we want the coords
    result = xr.DataArray(values, coords=coords)

    # if we stacked `da`, unstack result before returning
    return result.unstack(stack_dim) if stack_dim else result


def sparse_xr_dot(
    da1: xr.DataArray,
    da2: DataSetOrArray,
    debug: bool = False,
    broadcast_dims: Optional[Sequence[str]] = None,
) -> DataSetOrArray:
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
            if j in (i, 1):
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
