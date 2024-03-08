#!/usr/bin/env python
import functools
import json
import re
from pathlib import Path
from typing import Any, Literal, Optional, TypeVar

import pandas as pd
import xarray as xr
from xarray.core.common import DataWithCoords


# type for xr.Dataset *or* xr.DataArray
DataSetOrArray = TypeVar("DataSetOrArray", bound=DataWithCoords)


var_pat = re.compile(r"\s*[a-z]+ ([a-zA-Z_]+)\(.*\)")
attr_pat = re.compile(r"\s+([a-zA-Z_]+):([a-zA-Z_]+)\s*=\s*([^;]+)")


def get_data_var_attrs(template_file: str, species: Optional[str] = None) -> dict[str, dict[str, Any]]:
    """Extract data variable attributes from template file."""
    attr_dict: dict[str, Any] = {}

    with open(template_file, "r") as f:
        in_vars = False
        for line in f.readlines():
            if line.startswith("variables"):
                in_vars = True
            if in_vars:
                if m := var_pat.match(line):
                    attr_dict[m.group(1)] = {}
                if (m := attr_pat.match(line)) is not None and "FillValue" not in m.group(2):
                    val = m.group(3).strip().strip('"')

                    if species is not None:
                        val = val.replace("<species>", species)

                    attr_dict[m.group(1)][m.group(2)] = val

    return attr_dict


def make_global_attrs(
    output_type: Literal["flux", "conc"],
    author: str = "OpenGHG",
    species: str = "inert",
    domain: str = "EUROPE",
    apriori_description: str = "EDGAR 8.0",
    history: Optional[str] = None,
    comment: Optional[str] = None,
) -> dict[str, str]:
    global_attrs = {}
    global_attrs["title"] = (
        "Observed and simulated atmospheric concentrations"
        if output_type == "conc"
        else "Flux estimates: spatially-resolved and by country"
    )
    global_attrs.update(
        author=author,
        source="processed NAME(8.0) model output",
        transport_model="NAME",
        transport_model_version="NAME III (version 8.0)",
        met_model="UKV",
        species=species,
        domain=domain,
        inversion_method="RHIME",
        apriori_description=apriori_description,
        publication_acknowledgements="Please acknowledge ACRG, University of Bristol, in any publication that uses this data.",
    )
    global_attrs["history"] = history if history is not None else ""
    global_attrs["comment"] = comment if comment is not None else ""

    return global_attrs


def add_variable_attrs(
    ds: xr.Dataset, attrs: dict[str, dict[str, Any]], units: Optional[float] = None
) -> xr.Dataset:
    """Update data variables and coordinates of Dataset based on attributes dictionary.

    If `units` provided, data variables with "units" attribute will be rescaled by `units`. This is to convert e.g.
    from 1e-9 mol/mol to mol/mol.
    """
    for k, v in attrs.items():
        if k in ds.data_vars:
            if units is not None and "units" in v and "mol/mol" in v["units"]:
                ds[k] = units * ds[k]
            ds[k].attrs = v
        elif k in ds.coords:
            ds.coords[k].attrs = v

    return ds


@functools.lru_cache
def get_iso3166_codes() -> dict[str, Any]:
    """Load dictionary mapping alpha-2 country codes to other country information."""
    paris_formatting_path = Path(__file__).parent
    with open(paris_formatting_path / "iso3166.json", "r", encoding="utf8") as f:
        iso3166 = json.load(f)
    return iso3166


def get_country_code(
    x: str, iso3166: Optional[dict[str, dict[str, Any]]] = None, code: Literal["alpha2", "alpha3"] = "alpha3"
) -> str:
    """Get alpha-2 or alpha-3 (default) country code given the name of a country."""
    if iso3166 is None:
        iso3166 = get_iso3166_codes()

    # first try to match long names
    for v in iso3166.values():  # type: ignore
        if x.lower() == v["iso_long_name"].lower():
            return v[code]

    # next try to match unofficial names
    for v in iso3166.values():  # type: ignore
        if any(x.lower() == name.lower() for name in v["unofficial_names"]):
            return v[code]

    # next try to match substrings...
    for v in iso3166.values():
        names = [v["iso_long_name"].lower()] + [name.lower() for name in v["unofficial_names"]]
        if any(x.lower() in name for name in names):
            return v[code]

    # if no matches are found, return x
    return x


def convert_time_to_unix_epoch(x: DataSetOrArray, units: str = "1s") -> DataSetOrArray:
    """Convert `time` coordinate of xarray Dataset or DataArray to number of "units" since
    1 Jan 1970 (the "UNIX epoch").
    """
    return x.assign_coords(time=(pd.DatetimeIndex(x.time) - pd.Timestamp("1970-01-01")) // pd.Timedelta(units))
