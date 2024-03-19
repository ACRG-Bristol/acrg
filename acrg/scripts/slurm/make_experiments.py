import argparse
from pathlib import Path
from typing import Union

from helpers import flatten, make_dates_df, update_ini_file
from make_slurm_array import make_script

try:
    import tomllib
except ImportError:
    import pip._vendor.tomli as tomllib


def main(toml_path: Union[str, Path]) -> None:
    with open(toml_path, "rb") as f:
        conf = tomllib.load(f)

    setup = conf["setup"]

    kwargs_flat = flatten(conf["kwargs"])

    names_dict = {}
    for k, v in conf["kwargs"].items():
        if k in conf["names"]:
            names_dict[k] = conf["names"][k]
        elif not isinstance(v, list):
            names_dict[k] = None
        else:
            names_dict[k] = list(map(str, v))

    def make_name(x: dict) -> str:
        return "_".join([f"{k}_{v}" for k, v in x.items() if v is not None])

    names_flat = [make_name(x) for x in flatten(names_dict)]

    dates_df = make_dates_df(**conf["dates"])

    ini_template_path = Path(setup["ini_file"])
    job_root_path = Path(setup["job_root"])

    python_venv = setup.get("python_venv", None)
    conda_venv = setup.get("conda_venv", None)

    for i, (name, kwargs) in enumerate(zip(names_flat, kwargs_flat)):
        # set up directory to hold ini. slurm scipt, config, results, logs
        job_name = setup["job_name"]
        out_name = setup["out_prefix"] + "_" + name

        out_path = job_root_path / out_name

        if out_path.exists():
            j = 0
            while out_path.exists():
                new_out_name = out_name + str(j)
                out_path = job_root_path / new_out_name
                j += 1

        out_path.mkdir()

        # write config file with dates
        config_path = out_path / "inversion_dates.txt"
        dates_df.to_csv(config_path, sep="\t")

        # write ini file with updated kwargs
        ini_out_path = out_path / f"{job_name}.ini"
        updated_ini = update_ini_file(ini_template_path, new_kwargs=kwargs)

        with open(ini_out_path, "w") as f:
            f.writelines(updated_ini)

        # make slurm script
        slurm_script = make_script(
            job_name=job_name,
            config_file=config_path,
            ini_file=ini_out_path,
            out_dir=out_path,
            log_path=out_path,
            conda_env=conda_venv,
            python_venv=python_venv,
            n_array_jobs=conf["dates"]["n_periods"],
            n_kwargs=2,
            **conf["slurm"],
        )

        with open(out_path / "slurm.sh", "w") as f:
            f.write(slurm_script)

        with open(out_path / "readme.txt", "w") as f:
            f.write("Experiment using parameters:\n\n")
            for k, v in kwargs:
                f.write(f"{k}: {v}\n")

        print(f"Array config file and inversion wrapper for experiment {i} written to {out_path}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="experiment_maker",
        description="Make ini files and scripts for running several months of inversion on SLURM, possibly varying kwargs.",
    )
    parser.add_argument("toml_file")

    args = parser.parse_args()

    main(args.toml_file)
