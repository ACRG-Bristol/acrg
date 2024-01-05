#!/usr/bin/env python
"""
This script makes .ini files
"""
import argparse
from collections import defaultdict
import configparser
import csv
import json
from typing import Any, Optional


def get_paris_site_info(
    sites: Optional[list[str]] = None, filename: Optional[str] = None
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
            "TOB",
            "TOH",
            "TRN",
            "UTO",
            "WAO",
        ]

    if filename is None:
        filename = "site_params.csv"

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
    ini_template: Optional[str] = None,
    site_info: Optional[dict[str, list[Optional[str]]]] = None,
    output_name: Optional[str] = None,
) -> configparser.ConfigParser:
    """Load template ini file and update species and output paths, as well as sites,
    averaging_period, inlet, instrument, and fp_height with options provided by site_info
    """
    config = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))

    if ini_template is None:
        ini_template = "test.ini"

    config.read(ini_template)

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
    output_ini_file: str,
    species: str,
    root: str,
    ini_template: Optional[str] = None,
    site_params_file: Optional[str] = None,
    output_name: Optional[str] = None,
    sites: Optional[list[str]] = None,
) -> None:
    """Main script.

    Args:

    """
    if not output_ini_file.endswith(".ini"):
        output_ini_file = output_ini_file + ".ini"

    site_info = get_paris_site_info(sites, filename=site_params_file)

    kwargs = {}
    if ini_template is not None:
        kwargs["ini_template"] = ini_template
    if output_name is not None:
        kwargs["output_name"] = output_name

    config = make_ini_configparser(species, root, site_info=site_info, **kwargs)

    with open(output_ini_file, "w") as f:
        config.write(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("species", type=str)
    parser.add_argument("root", type=str)
    parser.add_argument("--sites", type=(lambda x: x.split(",")), default=None)
    parser.add_argument("--hbmcmc-output-name", type=str, default=None, help="name for hbmcmc output files")
    parser.add_argument(
        "--output-ini-file",
        type=str,
        default="nameless.ini",
        help="file name (with path) to save ini file output",
    )
    parser.add_argument("--ini-template", type=str, default=None, help="ini file to use as template")
    parser.add_argument(
        "--site-params",
        type=str,
        help="csv with site params; csv column names should be: sites, averaging_period, inlet, instrument, fp_height",
    )

    args = parser.parse_args()

    print("Command line args: ", vars(args))

    main(
        args.output_ini_file,
        args.species,
        args.root,
        sites=args.sites,
        output_name=args.hbmcmc_output_name,
        site_params_file=args.site_params,
        ini_template=args.ini_template,
    )
