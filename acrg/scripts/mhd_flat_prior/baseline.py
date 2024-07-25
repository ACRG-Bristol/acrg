#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Optional

import numpy as np
from openghg.types import SearchError
from openghg.standardise import standardise_bc
from openghg.retrieve import get_footprint, get_obs_surface
import pandas as pd
import xarray as xr


def create_dates(initial_year: int, n_years: int, freq: str = "MS") -> list[tuple[str, str]]:
    """Create iterator over tuples (start_date, end_date), where end_date - start_date is
    one month, and the start dates begin at the initial year, and continue for n_years at
    the given frequency.

    Use "QS" for quarterly frequency (e.g. for testing) and "MS" for (start of) monthly frequency.
    """
    start = pd.Timestamp(initial_year, 1, 1)
    dr = pd.date_range(start=start, end=(start + pd.DateOffset(years=n_years, months=-1)), freq=freq)

    sdr = [str(d).split("T")[0] for d in dr.values]
    dr_end = dr + pd.DateOffset(months=1)
    sdr_end = [str(d).split("T")[0] for d in dr_end.values]

    return list(zip(sdr, sdr_end))


def get_baseline(species: str, start_date: str, end_date: str, inlet: Optional[str] = None) -> dict[str, Any]:
    """Calculate baseline mole fraction using MHD observations when the wind direction is 180 to 300 degrees."""
    if inlet is None:
        obs = get_obs_surface(species=species, site="mhd", start_date=start_date, end_date=end_date).data
        wind = get_footprint(
            site="mhd", domain="europe", start_date=start_date, end_date=end_date
        ).data.wind_direction
    else:
        obs = get_obs_surface(
            species=species, site="mhd", start_date=start_date, end_date=end_date, inlet=inlet
        ).data
        wind = get_footprint(
            site="mhd",
            domain="europe",
            start_date=start_date,
            end_date=end_date,
            inlet=inlet,
        ).data.wind_direction

    wind = wind.reindex_like(obs, method="nearest")

    baseline = obs.where(wind < 300, drop=True).where(wind > 180, drop=True)

    result = dict(
        baseline=float(baseline.mf.mean().values),
        baseline_std=float(baseline.mf.std().values),
        percent_baseline=len(baseline.mf) / float(len(obs.mf)),
    )

    quantiles = np.quantile(obs.mf.values, np.linspace(0.1, 1.0, 9, endpoint=False))
    result.update({"q" + str(10 * (i + 1)): float(q) for i, q in enumerate(np.nditer(quantiles))})

    return result


def make_baseline_df(
    species: str, initial_year: int, n_years: int, inlet: Optional[str] = None, freq: Optional[str] = None
) -> pd.DataFrame:
    """Create DataFrame with MHD baseline mean and std over each month starting at 1 Jan `initial_year`,
    for `n_years` with frequency `freq` (default freq. is quarterly).

    You might need to specify the inlet height for some species (e.g. CH4); it is usually "10m", but might
    be "24m" (or something else).
    """
    dates = create_dates(initial_year, n_years, freq) if freq else create_dates(initial_year, n_years)

    results = []
    for start, end in dates:
        try:
            baseline = get_baseline(species, start, end, inlet)
        except (SearchError, AttributeError) as e:
            print(f"Error for start {start}: {e}")
            baseline = {"baseline": np.NaN, "baseline_std": np.NaN, "percent_baseline": np.NaN}
            baseline.update({"q" + str(10 * i): np.NaN for i in range(1, 10)})

        result = {"time": start}
        result.update(baseline)
        results.append(result)

    df = pd.DataFrame(results)
    df["time"] = pd.to_datetime(df["time"])
    df = df.set_index("time")

    return df


def create_flat_bc_prior(
    species: str,
    baseline_values: pd.Series,
    bc_template_file: Optional[str] = None,
    units: float = 1e-9,
    author: str = "aramsden",
    creation_method: Optional[str] = None,
) -> xr.Dataset:
    """Create a flat BC prior with the given baseline values based on a template file."""
    if bc_template_file is None:
        bc_template_file = "/group/chem/acrg/LPDM/bc/EUROPE/ch4_EUROPE_201901.nc"

    if creation_method is None:
        title = f"Uniform boundary conditions for {species}."
    else:
        title = f"Uniform boundary conditions for {species} created using {creation_method}."

    bc_template = xr.open_dataset(bc_template_file)
    bc_template = xr.ones_like(bc_template.drop_vars("time"))

    baseline_da = xr.DataArray(baseline_values)

    bc = bc_template * baseline_da * units
    bc.attrs["title"] = title
    bc.attrs["species"] = species
    bc.attrs["units"] = f"{units} mol / mol"
    bc.attrs["author"] = author
    bc.attrs["date_created"] = str(np.datetime64("now"))

    return bc


def main(
    species: str,
    initial_year: int,
    n_years: int,
    units: float,
    author: str = "OpenGHG cloud",
    creation_method: str = "MHD obs when wind direction between 180 and 300 degrees",
    freq: str = "MS",
    inlet: Optional[str] = None,
    output_dir: str = "/group/chem/acrg/LPDM/bc/EUROPE/paris_flat",
    uncert_output_dir: Optional[str] = None,
    standardise: bool = False,
    bc_input: Optional[str] = None,
    store: Optional[str] = None,
) -> None:
    baseline_df = make_baseline_df(species, initial_year, n_years, freq=freq, inlet=inlet).ffill()
    bc_ds = create_flat_bc_prior(
        species, baseline_df["baseline"], units=units, author=author, creation_method=creation_method
    )

    file_name = f"{species}_EUROPE_flat_{initial_year}-{initial_year + n_years}.nc"
    output_path = Path(output_dir)
    if not output_path.exists():
        output_path.mkdir(parents=True)
    file_path = output_path / file_name
    bc_ds.to_netcdf(file_path)

    if uncert_output_dir:
        uncert_path = Path(uncert_output_dir)
        if not uncert_path.exists():
            uncert_path.mkdir(parents=True)
        uncert_file_path = uncert_path / f"{species}_BC_uncert_{initial_year}-{initial_year + n_years}.csv"
        uncert = baseline_df["baseline_std"] / baseline_df["baseline"]
        uncert.to_csv(uncert_file_path)

    if standardise:
        bc_input = bc_input if bc_input is not None else "flat"
        standardise_bc(file_path, species, bc_input, "EUROPE", period=freq, store=store)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("species", type=str)
    parser.add_argument("initial_year", type=int)
    parser.add_argument("n_years", type=int)
    parser.add_argument("units", type=float)
    parser.add_argument("--author", type=str)
    parser.add_argument("--creation-method", type=str)
    parser.add_argument("--freq", type=str)
    parser.add_argument("--inlet", type=str)
    parser.add_argument("-o", "--output-dir", type=str)
    parser.add_argument("-u", "--uncert-output-dir", type=str)
    parser.add_argument("-s", "--standardise", default=False, action="store_true")
    parser.add_argument("--bc-input", type=str)
    parser.add_argument("--store", type=str)

    args = vars(parser.parse_args())

    main(**{k: v for k, v in args.items() if v is not None})
