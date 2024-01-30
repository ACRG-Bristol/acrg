#!/usr/bin/env python3
"""
Script to create array job scripts to run the same inversion over several months.

NOTE: You need to set `outputpath` in your .ini file to whatever you specify in the
script. This can't be modified through the command line with the current version of
`run_hbmcmc.py`
"""
import argparse
import csv
import json
from pathlib import Path
import textwrap
from typing import Optional, Union

import pandas as pd


# where should jobs be stored?
DEFAULT_JOB_ROOT = Path("/user/work/bm13805")
DEFAULT_INVERSIONS_PATH = Path("/user/work/bm13805/openghg_inversions")


def create_array_config_df(year: int, n_months: int, initial_month: int = 1, **kwargs) -> pd.DataFrame:
    start = pd.Timestamp(year, initial_month, 1)
    dr = pd.date_range(start=start, end=(start + pd.DateOffset(months=n_months)), freq="M")

    # dr contains last days of each month; we want the end dates to be the first day of the next month
    end_dates = dr + pd.DateOffset(days=1)
    start_dates = end_dates - pd.DateOffset(months=1)

    array_task_id = list(range(1, n_months + 1))
    result = pd.DataFrame(
        data={"array_task_id": array_task_id, "start_date": start_dates, "end_date": end_dates}
    ).set_index("array_task_id")

    if kwargs:
        for k, v in kwargs.items():
            result[k] = [v] * n_months
            result[k + "_code"] = [list(range(1, len(v) + 1))] * n_months

        kw_cols = list(kwargs.keys())
        kw_encoding = [k + "_code" for k in kw_cols]

        # expand parameter lists to get all keyword combinations
        for k in kw_cols:
            result = result.explode([k, k + "_code"], ignore_index=True).rename_axis(result.index.name)

        # combine kw cols into strings that will be parsed by json.loads to a dictionary
        result["kwargs"] = result[kw_cols].apply(
            lambda x: json.dumps({k: x[i] for i, k in enumerate(kw_cols)}), axis=1
        )
        result["kwargs_dir_name"] = result[kw_encoding].apply(
            lambda x: "_".join([k + "_" + str(x.loc[kw_encoding[i]]) for i, k in enumerate(kw_cols)]), axis=1
        )
        result = result.drop(columns=(kw_cols + kw_encoding))

    result.index = result.index + 1  # slurm array jobs start at 1

    return result


def make_script(
    job_name: str,
    config_file: Path,
    ini_file: Path,
    out_dir: Path = DEFAULT_JOB_ROOT,
    log_path: Path = DEFAULT_JOB_ROOT,
    n_cpu: int = 4,
    time: str = "02:00:00",
    mem: str = "20gb",
    conda_env: Optional[str] = None,
    python_venv: Optional[str] = None,
    n_months: Optional[int] = None,
) -> str:
    """Make script to run array job on SLURM.

    TODO: need to set output paths for .out, .err, and .ini
    """
    # Get number of jobs in array
    with open(config_file, "r") as f:
        keywords = f.readline()[:-1].split("\t")  # remove trailing \n and split
        n_array_jobs = len(f.readlines())

    # Make strings for each section of wrapper script
    header_str = f"""\
    #!/bin/bash
    #SBATCH --job-name={job_name}
    #SBATCH --output={log_path / job_name}_%a.out
    #SBATCH --error={log_path / job_name}_%a.err
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=1
    #SBATCH --cpus-per-task={n_cpu}
    #SBATCH --time={time}
    #SBATCH --mem={mem}
    #SBATCH --account=chem007981
    #SBATCH --array=1-{n_array_jobs}
    """

    if conda_env:
        env_str = f"""\
        # Set up Python environment
        module purge
        module load lang/python/anaconda
        module load tools/git
        eval "$(conda shell.bash hook)"
        conda activate {conda_env}
        """
        env_log_str = f"""\
        if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
            mkdir -p {log_path}
            conda list >> {log_path / "packages.txt"}
        fi
        """
    elif python_venv:
        env_str = f"""\
        # Set up Python environment
        module purge
        module load lang/python/anaconda
        module load tools/git
        source {python_venv}/bin/activate
        """
        env_log_str = f"""\
        if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
            mkdir -p {log_path}
            pip list >> {log_path / "packages.txt"}
        fi
        """
    else:
        raise ValueError("One of `conda_env` or `python_venv` must be supplied.")

    # check numpy config and see if pymc loads without warnings
    # TODO

    # record what branch and commit was used for these inversions
    branch_str = f"""\
    inversions_path=$(pip list | grep "openghg-inversion" | awk '{{ print $3 }}')
    git_branch=$(git -C $inversions_path status | awk 'NR==1{{ print $3 }}')
    git_commit=$(git -C $inversions_path log --oneline -n 1 $git_branch | awk '{{ print $1 }}')
    echo "Using commit $git_commit on branch $git_branch in repo $inversions_path" >> {log_path / "git_info.txt"}
    """

    # set up parameters to be passed to run_hbmcmc
    param_str = f"""\
    # Specify the path to the config file
    config={config_file}
    nkwargs={len(keywords)}

    # Extract start and end dates for the current $SLURM_ARRAY_TASK_ID
    start=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\\t' '$1==array_task_id {{print $2}}' $config)
    end=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\\t' '$1==array_task_id {{print $3}}' $config)

    # Keyword arg handling
    if [ $nkwargs -gt 3 ]; then
        # Extract keyword args for the current $SLURM_ARRAY_TASK_ID
        kwargs=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\\t' '$1==array_task_id {{print $4}}' $config)

        # Extract output dir name corresponding to keyword args
        child_out_dir=$(awk -v array_task_id=$SLURM_ARRAY_TASK_ID -F'\\t' '$1==array_task_id {{print $5}}' $config)
    fi
    """

    # run the inversion
    command_str = f"""\

    # Print a message with current $SLURM_ARRAY_TASK_ID, start date, and end date
    echo "Running array task ${{SLURM_ARRAY_TASK_ID}}: start date ${{start}}, end date is ${{end}}."

    # Check if there are kwargs
    if [ $nkwargs -le 3 ]; then
        python ${{inversions_path}}/openghg_inversions/hbmcmc/run_hbmcmc.py "${{start}}" "${{end}}" -c {ini_file} --output-path={out_dir}
    else
        echo "Keyword args: ${{kwargs}}"
        echo "Output directory name: ${{child_out_dir}}"
        python ${{inversions_path}}/openghg_inversions/hbmcmc/run_hbmcmc.py "${{start}}" "${{end}}" -c {ini_file} --kwargs="${{kwargs}}" --output-path="{out_dir}/${{child_out_dir}}"
    fi
    """

    # remove leading whitespace and combine sections
    sections = [header_str, env_str, env_log_str, branch_str, param_str, command_str]
    script_str = "\n".join([textwrap.dedent(section) for section in sections])
    return script_str


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="array_job_maker",
        description="Make script for running several months of inversion on SLURM.",
    )
    parser.add_argument("-J", "--jobname")
    parser.add_argument(
        "--ini-file",
        help=".ini file containing RHIME hbmcmc config.",
    )
    parser.add_argument("--ini-path", help="path to .ini file", default=DEFAULT_JOB_ROOT)
    parser.add_argument(
        "-r",
        "--job-root",
        help="directory JOB_ROOT where directories for individual jobs will be stored.",
        type=Path,
        default=DEFAULT_JOB_ROOT,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="directory within JOB_ROOT to store files related to this job.",
        default=None,
    )
    parser.add_argument("--conda-env", help="name of conda environment to use", default=None)
    parser.add_argument("--python-venv", help="path to python virtual environment to use", default=None)

    # dates for inversion
    parser.add_argument("-y", "--year", type=int)
    parser.add_argument("-n", "--nmonths", type=int)
    parser.add_argument("--initial-month", type=int, default=1)

    # keyword arguments
    parser.add_argument(
        "--kwargs", type=json.loads, help="additional keyword arguments to pass to hbmcmc, in json format."
    )

    # params for compute requirements
    parser.add_argument("--ncpu", type=int, help="number of CPUs per job", default=4)
    parser.add_argument("--mem", help="amount of memory for jobs", default="20gb")
    parser.add_argument("--time", help="length of inversion runs", default="02:00:00")

    args = vars(parser.parse_args())  # get args in dictionary
    print(args)
    job_name = args["jobname"]
    out_path = args["job_root"] / job_name.replace(" ", "_")
    if args["outdir"]:
        out_path = out_path / args["outdir"]
    out_path.mkdir(parents=True, exist_ok=True)

    # create config TSV
    kwargs = args["kwargs"] if args["kwargs"] else {}
    array_config_df = create_array_config_df(
        args["year"],
        args["nmonths"],
        args["initial_month"],
        **kwargs,
    )
    config_path = out_path / "config.txt"
    array_config_df.to_csv(config_path, sep="\t", doublequote=False, quotechar="'", quoting=csv.QUOTE_NONE)

    # make wrapper script
    wrapper_path = out_path / "inversion_wrapper.sh"
    script_str = make_script(
        job_name=job_name,
        config_file=config_path,
        ini_file=Path(args["ini_path"]) / args["ini_file"],
        n_cpu=args["ncpu"],
        mem=args["mem"],
        time=args["time"],
        conda_env=args["conda_env"],
        python_venv=args["python_venv"],
        out_dir=out_path,
        log_path=(out_path / "logs"),
    )
    with open(wrapper_path, "w") as f:
        f.write(script_str)

    print(f"Array config file and inversion wrapper written to {out_path}.")
