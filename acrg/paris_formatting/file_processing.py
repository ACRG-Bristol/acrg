#!/usr/bin/env python
from pathlib import Path
from typing import Optional, Union

import numpy as np
import xarray as xr


def get_netcdf_files(dir: Union[str, Path], filename_search: Optional[str] = None) -> list[Path]:
    """Get list of paths to netCDF files, optionally filtered by `filename_search` string."""
    if isinstance(dir, str):
        dir = Path(dir)

    if filename_search is None:
        file_list = sorted(dir.glob("*.nc"))
    else:
        file_list = sorted(dir.glob(f"*{filename_search}*.nc"))

    return file_list


def get_rhime_outs(files: list[Path]) -> list[xr.Dataset]:
    """Load RHIME output files as xr.Datasets and process some variables.

    NOTE: some coordinates in the RHIME output files seem to be duplicated (e.g. paramnum and nlatent?)
    TODO: Check/fix this in inversions.
    """
    rename_dict = dict(stepnum="draw", paramnum="nlatent", numBC="nBC", measurenum="nmeasure")
    outs = [
        xr.open_dataset(file)
        .rename_vars(rename_dict)
        .drop_dims(["nUI", "nlatent"])
        .swap_dims({"steps": "draw"})
        for file in files
    ]

    for ds in outs:
        if len(ds.xtrace.coords) < 2:
            ds["xtrace"] = ds.xtrace.assign_coords(nparam=("nparam", np.arange(ds.dims["nparam"])))

    return outs


def nmeasure_to_site_and_time(ds: xr.Dataset) -> xr.Dataset:
    """Rough function extracted from flux_output_format.py"""
    time_vars = [dv for dv in ds.data_vars if "nmeasure" in list(ds[dv].coords)]
    time_vars = [
        tv
        for tv in time_vars
        if all(x not in str(tv) for x in ["mode", "median", "off", "sensitivity", "sigma"])
    ]

    ds = ds[time_vars].drop_dims("nUI")

    site_nums = np.arange(ds.siteindicator.max() + 1)

    result = (
        xr.concat(
            [
                ds.where(ds.siteindicator == site_num, drop=True)
                .expand_dims({"nsite": [site_num]})
                .assign_coords(nmeasure=ds.Ytime.where(ds.siteindicator == site_num, drop=True))
                .rename_vars(nmeasure="time")
                for site_num in site_nums
            ],
            dim="nsite",
        )
        .swap_dims(nmeasure="time")
        .drop_vars(["Ytime", "siteindicator"])
    )
    return result
