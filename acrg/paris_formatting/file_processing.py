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
    rename_dict = dict(stepnum="steps", paramnum="nlatent", numBC="nBC", measurenum="nmeasure", UInum="nUI")
    outs = [xr.open_dataset(file).rename_vars(rename_dict) for file in files]

    for ds in outs:
        if len(ds.xtrace.coords) < 2:
            ds["xtrace"] = ds.xtrace.assign_coords(nparam=("nparam", np.arange(ds.dims["nparam"])))

    return outs
