#!/usr/bin/env python
"""
This script searches for observation and footprint data, as specified by a .ini file.

It will show you the search results as well as the number of days of missing data for
each result.

Example usage:

python scenario_data_checker.py "/user/work/bm13805/ini_files/test.ini" "2018-01-01" "2022-12-31"

You need to be in a virtual environment with openghg and openghg_inversions installed.

The script prints the results by default, with missing data presented by year.
The `-v` flag suppresses these print statements.

The `-s` flag will save the outputs to `<path>/<ini file name>_search_results.csv` and
`<path>/<ini file name>_missing_data.csv`, where `<path>/<ini file name>.ini` is the first
argument to the script.

You can aggregate missing days over other frequencies. For monthly, pass the argument `--freq "M"`.
"""
from argparse import ArgumentParser
import json
from pathlib import Path
from typing import Optional

from openghg.objectstore import get_readable_buckets
from openghg.objectstore.metastore import open_metastore
from openghg.retrieve import search_surface, search_footprints, search_flux, search_bc
from openghg.util import clean_string, format_inlet
import openghg_inversions.config.config as oiconfig
import pandas as pd


def get_search_result_dateranges(res: pd.Series) -> list[str]:
    """Get the daterange strings from the datasource associated with a row
    from `openghg.retrieve.search().results`
    """
    store_path = Path(res["object_store"])
    datasource_path = store_path / "datasource" / "uuid" / (res["uuid"] + "._data")

    with open(datasource_path, "r") as f:
        metadata = json.load(f)

    return list(metadata["data_keys"]["latest"]["keys"])


def pd_daterange_from_str(
    dr: str, start: Optional[str] = None, end: Optional[str] = None
) -> pd.DatetimeIndex:
    _start, _end = dr.split("_")
    _start = _start[:10]  # extract year-month-day
    _end = _end[:10]  # extract year-month-day

    if start is not None and start > _start:
        _start = start

    if end is not None and end < _end:
        _end = end

    return pd.date_range(_start, _end)


def get_search_result_series(
    res: pd.Series, start: Optional[str] = None, end: Optional[str] = None
) -> pd.Series:
    """Return 0-1 valued series covering start to end (or largest range covered by
    search results).

    Use like:
    `search_results.apply(lambda x: get_search_results_series(x, "2018-01-01", "2022-12-31"), axis=1)`
    """
    start = start if start is not None else res["start_date"][:10]
    end = end if end is not None else res["end_date"][:10]
    index = pd.date_range(start, end)
    result = pd.Series(0, index=index)

    drs = get_search_result_dateranges(res)

    for dr in drs:
        pd_dr = pd_daterange_from_str(dr, start, end)
        result[pd_dr] = 1

    return result


def search_by_ini(ini_path: str) -> pd.DataFrame:
    """Search for surface and footprint data as specified by an .ini file.

    This searches for all surface and footprint data that would be used in an inversion.
    """
    param = dict(oiconfig.extract_params(ini_path))
    rbuckets = get_readable_buckets()
    search_results = []

    tmp = param["obs_store"]
    obs_stores = tmp if isinstance(tmp, list) else [tmp]

    tmp = param["footprint_store"]
    footprint_stores = tmp if isinstance(tmp, list) else [tmp]

    # get surface results
    species = clean_string(param["species"])

    for store in obs_stores:
        try:
            surface_bucket = rbuckets[store]
        except KeyError:
            raise ValueError(f"Store {store} not found. Run `openghg --quickstart` to add store.")

        with open_metastore(bucket=surface_bucket, data_type="surface", mode="r") as metastore:
            for site, inlet, instrument in zip(param["sites"], param["inlet"], param["instrument"]):
                sterms = dict(species=species, site=clean_string(site))

                if inlet is not None:
                    sterms["inlet"] = format_inlet(inlet)
                if instrument is not None:
                    sterms["instrument"] = clean_string(instrument)

                if result := metastore.search(sterms):
                    for r in result:
                        r["object_store"] = surface_bucket
                    search_results.extend(result)

    # get footprint results
    for store in footprint_stores:
        try:
            footprint_bucket = rbuckets[store]
        except KeyError:
            raise ValueError(f"Store {store} not found. Run `openghg --quickstart` to add store.")

        with open_metastore(bucket=footprint_bucket, data_type="footprints", mode="r") as metastore:
            for site, inlet in zip(param["sites"], param["fp_height"]):
                sterms = dict(site=clean_string(site), domain="europe")

                if inlet is not None:
                    sterms["inlet"] = format_inlet(inlet)

                if result := metastore.search(sterms):
                    for r in result:
                        r["object_store"] = footprint_bucket
                    search_results.extend(result)

    df = (
        pd.DataFrame.from_dict(search_results)
        .loc[
            :, ["site", "data_type", "inlet", "instrument", "start_date", "end_date", "uuid", "object_store"]
        ]
        .sort_values("site")
        .rename_axis(index="id")
        .set_index(["site", "data_type"], append=True)
        .reorder_levels([1, 2, 0], axis=0)
    )

    return df


def get_data_times_df(search_results: pd.DataFrame, start: str, end: str) -> pd.DataFrame:
    """Get data frame whose columns are the index of `search_results` with daily date range index
    from start to end. The entries are 1 if data is present, 0 if not.

    Args:
        search_results: OpenGHG `SearchResults.results` or the result of `search_by_ini`.
    """
    result = search_results.apply(lambda x: get_search_result_series(x, start, end), axis=1).T.rename_axis(
        index="date"
    )
    return result


def get_missing_data_times_df(
    data_times: pd.DataFrame, resample_to: str = "Y", group_results=False
) -> pd.DataFrame:
    """Report total days of missing data over time periods specified by `resample_to`.

    Args:
        data_times: result of `get_data_times_df`
    """
    if group_results:
        result = ((1 - data_times).T.groupby(["site", "data_type"]).sum() > 0).T.resample(resample_to).sum().T
    else:
        result = (1 - data_times).T.resample(resample_to).sum().T
    return result


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("ini_file", type=str)
    parser.add_argument("start", type=str)
    parser.add_argument("end", type=str)
    parser.add_argument("--freq", type=str, default="Y")
    parser.add_argument(
        "-s",
        "--save-output",
        default=False,
        action="store_true",
        help="file name to use for output; don't give a file extension.",
    )
    parser.add_argument(
        "-v", "--silent", default=False, action="store_true", help="suppress print statements"
    )

    args = parser.parse_args()

    search_results = search_by_ini(args.ini_file)
    df = get_data_times_df(search_results, args.start, args.end)
    missing_df = get_missing_data_times_df(df, resample_to=args.freq, group_results=True)
    search_results["total_missing"] = get_missing_data_times_df(df, resample_to=args.freq).sum(axis=1)

    if not args.silent:
        with pd.option_context("display.max_rows", None, "display.max_columns", None):
            print("\nSearch results:\n:", search_results[["inlet", "instrument", "uuid", "total_missing"]])
            print(f"Missing data by {args.freq}:\n", missing_df)

    if args.save_output:
        search_results.to_csv((args.ini_file[:-4] + "_search_results.csv"))
        df.to_csv((args.ini_file[:-4] + "_data_times.csv"))
        missing_df.to_csv((args.ini_file[:-4] + "_missing_data.csv"))
