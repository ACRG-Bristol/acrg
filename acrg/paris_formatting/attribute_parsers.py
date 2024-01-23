#!/usr/bin/env python
import json
from pathlib import Path
import re
from typing import Any, Literal, Optional, Union

import xarray as xr

concentrations_template = "netcdf_template_concentrations_bm_edits.txt"
emissions_template = "netcdf_template_emissions_bm_edits.txt"

var_pat = re.compile(r"\s*[a-z]+ ([a-zA-Z]+)\(.*\)")
attr_pat = re.compile(r"\s+([a-zA-Z]+):([a-zA-Z_]+)\s*=\s*([^;]+)")


def get_data_var_attrs(template_file: str, drop_fill_value=True) -> dict[str, dict[str, Any]]:
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
                    attr_dict[m.group(1)][m.group(2)] = m.group(3).strip().strip('"')

    return attr_dict


def write_data_var_attrs(attr_dict: dict[str, Any], output_file: str) -> None:
    if not output_file.endswith(".json"):
        output_file = output_file + ".json"

    with open(output_file, "w") as f:
        json.dump(attr_dict, f)


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
    """Update data variables and coordinates of Dataset based on attributes dicitonary.

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
