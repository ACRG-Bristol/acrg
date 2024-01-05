#!/usr/bin/env python
import argparse
from collections import defaultdict
import configparser
import csv
import json
from typing import Any, Optional


def get_paris_params(filename: str = "model_scenario.params.json") -> dict[str, Any]:
    """Get params frozen from sample PARIS .ini file."""
    with open(filename, "r") as f:
        params = json.load(f)
    return params


def get_paris_site_info(
    sites: Optional[list[str]] = None, filename: str = "site_params.csv"
) -> dict[str, list[Optional[str]]]:
    """Get site code, averaging period, inlet, instrument, and fp height for selected sites."""
    if sites is None:
        sites = [
            "BIR",
            "BSD",
            "CBW",
            "CMN",
            "GAT",
            "HEI",
            "HEL",
            "HFD",
            "HPB",
            "HTM",
            "HUN",
            "JFJ",
            "KIT",
            "KRE",
            "LIN",
            "LUT",
            "MHD",
            "NOR",
            "OPE",
            "OXK",
            "PAL",
            "RGL",
            "SAC",
            "STE",
            "TAC",
            "TOH",
            "TRN",
            "UTO",
            "WAO",
        ]

    result = defaultdict(list)

    with open(filename, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if row[0] in sites:
                for k, v in zip(header, row):
                    result[k].append(v) if v else result[k].append(None)
    return result


def make_ini_configparser(
    species: str,
    root: str,
    site_info: Optional[dict[str, list[Optional[str]]]],
    output_name: Optional[str] = None,
) -> configparser.ConfigParser:
    config = configparser.ConfigParser()
    config.read("test.ini")

    config["INPUT.MEASUREMENTS"]["species"] = species
    config["INPUT.MEASUREMENTS"]["merged_data_dir"] = f"{root}/merged_data"

    if site_info:
        for k, v in site_info.items():
            if k == "fp_height":
                config["INPUT.PRIORS"][k] = str(v)
            else:
                config["INPUT.MEASUREMENTS"][k] = str(v)

    config["MCMC.OUTPUT"]["outputpath"] = f"{root}/{species}/hbmcmc_input_output"
    config["MCMC.OUTPUT"]["outputname"] = output_name if output_name else "hbmcmc_output"

    return config


def main(
    ini_file_name: str,
    species: str,
    root: str,
    output_name: Optional[str] = None,
    sites: Optional[list[str]] = None,
) -> None:
    if not ini_file_name.endswith(".ini"):
        output_name = output_name + ".ini"

    site_info = get_paris_site_info(sites)

    config = make_ini_configparser(species, root, site_info, output_name)

    with open(ini_file_name, "w") as f:
        config.write(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("species", type=str)
    parser.add_argument("root", type=str)
    parser.add_argument("--sites", type=(lambda x: x.split(",")), default=None)
    parser.add_argument("--output-name", type=str, default=None)
    parser.add_argument("-o", "--ini-file-name", type=str, default="nameless.ini")

    args = parser.parse_args()

    main(args.ini_file_name, args.species, args.root, sites=args.sites, output_name=args.output_name)
