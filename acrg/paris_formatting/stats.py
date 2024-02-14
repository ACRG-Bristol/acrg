from typing import Optional, Sequence

import numpy as np
import scipy
import xarray as xr


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

        if np.nanmax(row) > np.nanmin(row):
            xvals = np.linspace(np.nanmin(row), np.nanmax(row), 200)
            kde = scipy.stats.gaussian_kde(row).evaluate(xvals)
            return xvals[kde.argmax()]

        return np.nanmean(row)

    def func(arr):
        return np.apply_along_axis(func1d=mode_of_row, axis=-1, arr=arr)

    return xr.apply_ufunc(func, da, input_core_dims=[[sample_dim]], dask="parallelized")


def calculate_stats(
    ds: xr.Dataset,
    name: str,
    chunk_dim: str,
    chunk_size: int = 10,
    var_names: Optional[list[str]] = None,
    report_mode: bool = False,
    add_bc_suffix: bool = False,
) -> list[xr.Dataset]:
    output = []
    if var_names is None:
        var_names = list(ds.data_vars)
    for var_name in var_names:
        suffix = "apost" if "posterior" in var_name else "apriori"
        if add_bc_suffix:
            suffix += "BC"

        if report_mode:
            stats = [
                calc_mode(ds[var_name].dropna(dim="draw").chunk({chunk_dim: chunk_size}), sample_dim="draw")
                .compute()
                .rename(f"{name}{suffix}"),
            ]
        else:
            stats = [
                ds[var_name].mean("draw").rename(f"{name}{suffix}"),
                calc_mode(ds[var_name].dropna(dim="draw").chunk({chunk_dim: chunk_size}), sample_dim="draw")
                .compute()
                .rename(f"{name}{suffix}_mode"),
            ]
        stats.append(
            make_quantiles(ds[var_name].dropna(dim="draw"), sample_dim="draw").rename(f"q{name}{suffix}")
        )

        output.extend(stats)
    return output
