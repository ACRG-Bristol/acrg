import json
import re
from itertools import product
from pathlib import Path
from typing import Iterable, Literal, Optional, Union

import pandas as pd


def flatten(x: dict) -> list[dict]:
    """Flatten any iterable values in dictionary.

    For example:

    >>>flatten({"a": (1, 2), "b": "XY"})
    [{"a": 1, "b": "X"}, {"a": 1, "b": "Y"}, {"a": 2, "b": "X"}, {"a": 2, "b": "Y"}]

    >>>flatten({"a": (1, 2), "b": 3})
    [{"a": 1, "b": 3}, {"a": 2, "b": 3}]

    Args:
        x: dictionary to flatten

    Returns:
        list of dictionaries with same keys as input dictionary and no iterable values;
    all combinations of iterable values in the input dictionary appear as a value in
    some dictionary in the output.
    """
    keys, vals = zip(*x.items())  # get keys, values as tuples

    # all args of `product` must be Iterable, so convert non-Iterable values to lists of length 1
    def make_iterable(x):
        if isinstance(x, str) or not isinstance(x, Iterable):
            return [x]
        return x

    vals = map(make_iterable, vals)
    vals_prod = product(*vals)  # create tuples containing all combinations of values for each key

    return [dict(zip(keys, vp)) for vp in vals_prod]


def make_dates_df(
    year: int,
    n_periods: int,
    frequency: Literal["annual", "monthly"] = "annual",
    initial_month: int = 1,
    array_job_id: bool = True,
) -> pd.DataFrame:
    """Create a DataFrame containing `n_periods` start and end dates starting at the given
    `year` and initial month (`initial_month` defaults to 1).

    Args:
        year: year to start first period
        n_periods: number of periods (months or years) in result
        frequency: length of periods, "annual" or "monthly"
        initial_month: month to start first period (default 1)

    Returns:
        DataFrame containing columns for start and end dates.
    """
    if frequency == "annual":
        freq = "YS"
        n_years, n_months = n_periods, 0
        offset = pd.DateOffset(years=1)  # offset for start vs. end dates
    elif frequency == "monthly":
        freq = "MS"
        n_years, n_months = 0, n_periods
        offset = pd.DateOffset(months=1)  # offset for start vs. end dates
    else:
        raise ValueError(f"Frequency {frequency} not accepted.")

    start = pd.Timestamp(year, initial_month, 1)
    end = start + pd.DateOffset(years=n_years, months=n_months)  # type: ignore

    start_dates = pd.date_range(start, end, inclusive="left", freq=freq)
    end_dates = start_dates + offset

    dates_df = pd.DataFrame({"start_date": start_dates, "end_date": end_dates})

    if array_job_id:
        dates_df.index +=1
        dates_df = dates_df.rename_axis("array_job_id")

    return dates_df


def update_ini_file(ini_file: Union[str, Path], new_kwargs: Optional[dict] = None) -> list[str]:
    """Read ini file and update key-value pairs according to `new_kwargs`.

    If a key in `new_kwargs` matches a key in the .ini file, the value of that
    key is set to the value in `new_kwargs`.

    Any items in `new_kwargs` that are not matched are added to the end of the
    the result in a new section.

    Args:
        ini_file: path to .ini file
        new_kwargs: dictionary of new key-value pairs to use.

    Returns:
        list of lines of the updated ini file.
    """

    def make_kv_string(k, v):
        if isinstance(v, (dict, list, tuple)):
            return f"{k} = {json.dumps(v)}\n"
        elif isinstance(v, str):
            return f'{k} = "{v}"\n'
        return f"{k} = {v}\n"

    kv_pat = re.compile(r'([;#\s])*(\w+)\s*=\s*([-\[\] "/,\w]+)')
    kwargs_copy = new_kwargs.copy() if new_kwargs else {}
    conf_lines = []

    with open(ini_file, "r") as f:
        for line in f:
            if m := kv_pat.match(line):
                if (k := m.group(2)) in kwargs_copy:
                    v = kwargs_copy.pop(k)
                    new_line = make_kv_string(k, v)
                    conf_lines.append(new_line)
                else:
                    conf_lines.append(line)
            else:
                conf_lines.append(line)

    if kwargs_copy:
        conf_lines.extend(["\n", "[EXTRA.KWARGS]", "\n"])
        for k, v in kwargs_copy.items():
            conf_lines.append(make_kv_string(k, v))

    return conf_lines
