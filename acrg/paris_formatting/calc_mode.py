#!/usr/bin/env python3

import numpy as np
import openghg_inversions as oi
from openghg_inversions import utils
import pandas as pd
from scipy import stats
import sparse
from sparse import COO
import xarray as xr


def calc_mode(da: xr.DataArray, sample_dim="steps") -> xr.DataArray:
    """Calculate the (KDE smoothed) mode of a data array containing MCMC
    samples.
    """

    def mode_of_row(row):
        if np.nanmax(row) > np.nanmin(row):
            xvals = np.linspace(np.nanmin(row), np.nanmax(row), 200)
            kde = stats.gaussian_kde(row).evaluate(xvals)
            return xvals[kde.argmax()]
        else:
            return np.nanmean(row)

    def func(arr):
        return np.apply_along_axis(func1d=mode_of_row, axis=-1, arr=arr)

    result = xr.apply_ufunc(func, da, input_core_dims=[[sample_dim]])
