from __future__ import annotations
from dataclasses import dataclass, field
from functools import reduce
from operator import attrgetter, itemgetter
from pathlib import Path
from pprint import pprint
import re
from typing import Any, Callable, Iterator, Literal, Optional, TypeVar, Union

import pandas as pd

from helpers import flatten, make_dates_df, make_iterable, update_ini_file
from make_slurm_array import make_script

try:
    import tomllib
except ImportError:
    import pip._vendor.tomli as tomllib


# flag value for parameterising values by dates dervied from array job pd.DataFrame,
# which has "start_date" and "end_date" columns
skip_flags = ["array_row"]


def load_conf(toml_path: Union[str, Path]) -> dict:
    with open(toml_path, "rb") as f:
        conf = tomllib.load(f)
    return conf


@dataclass
class Dates:
    year: int
    n_periods: int
    frequency: Literal["monthly", "annual", "yearly"] = "annual"
    initial_month: int = 1
    array_job_id: bool = True

    def to_df(self) -> pd.DataFrame:
        return make_dates_df(self.year, self.n_periods, self.frequency, self.initial_month, self.array_job_id)


# regexes for toml keys and substitution parameters in config files
# The idea is that "<general.species>" will be replaced by the value of "species" in
# the "general" section of the config file.
#
# An exception: e.g. "<array_row.start_date.year>" is treated specially, since array jobs may span many years
# so this info is put into the slurm config.txt file used by the slurm script

# toml bare keys:
bare_toml_key_str = r"[A-Za-z0-9_-]+"
bare_toml_key_pat = re.compile(bare_toml_key_str)

# the following pattern captures toml keys: bare key(s), separated by "." + whitespace
# NOTE: it only captures the first and last groups, so once the full match is found,
# the bare keys can be extracted from the full key by the bare_toml_key_pat
toml_key_str = "(" + bare_toml_key_str + ")" + rf"(\s*\.\s*{bare_toml_key_str})*"
toml_key_pat = re.compile(toml_key_str)


def find_bare_keys(s: str) -> list[str]:
    """If s is a match to toml_key_pat, say "<general.species>", then
    this will return the list ["general", "species"]
    """
    return [m[0] for m in bare_toml_key_pat.finditer(s)]


sub_pat = re.compile(rf"<({toml_key_str})>")


@dataclass
class Substitution:
    start: int
    end: int
    match: str

    @property
    def keys(self) -> list[str]:
        return find_bare_keys(self.match)

    @property
    def key(self) -> str:
        return ".".join(find_bare_keys(self.match))

    @property
    def root(self) -> str:
        return self.keys[0]

    def apply(self, value: str, style: Literal["config", "ini"] | str = "config") -> str:
        # left = value[: self.start]
        # right = value[self.end :]
        # center = "{" + f"config.{self.root}" + "".join([f"[{key}]" for key in self.keys[1:]]) + "}"
        # return left + center + right
        old = value[self.start:self.end]
        if style in ["config", "ini"]:
            new = "{" + f"{style}.{self.root}" + "".join([f"[{key}]" for key in self.keys[1:]]) + "}"
        else:
            new = "{" + style + "}"
        return value.replace(old, new)

    def __hash__(self) :
        return hash(tuple(self.keys))


@dataclass
class Param:
    value: Any
    key: Optional[str] = None

    def __post_init__(self) -> None:
        self.vtype = type(self.value)

    def __str__(self) -> str:
        return str(self.value)

    def find_subsitutions(self):
        if not isinstance(self.value, str):
            return None

        subs = [Substitution(m.start(), m.end(), m[0]) for m in sub_pat.finditer(self.value)]
        return subs

    def _find_next_sub(self, skip: Optional[list[str]] = None) -> Optional[Substitution]:
        """Find next substitution whose root is not in `skip`."""
        subs = self.find_subsitutions()

        if subs is None:
            return None

        if skip is None:
            skip = []

        subs = [s for s in subs if not s.root in skip]

        if not subs:
            return None

        return subs[0]

    def _format(self, config: Config, skip: list[str] | None) -> None:
        if isinstance(self.value, str):
            # iteratively replace <key1.key2> with {config.key1[key2]} to set up for formatting below
            while s := self._find_next_sub(skip=skip):
                self.value = s.apply(self.value)

            self.value = self.value.format(config=config)

    def format(self, config: Config, skip: list[str] | None = None) -> None:
        """Recursively apply _format"""
        # print("formatting", self.key, self.value)
        if isinstance(self.value, list):
            tmp = [Param(x) for x in self.value]
            for x in tmp:
                x.format(config)
            self.value = [x.value for x in tmp]

        if isinstance(self.value, dict):
            tmp = {k: Param(v) for k, v in self.value.items()}
            for x in tmp.values():
                x.format(config)
            self.value = {k: v.value for k, v in tmp.items()}

        # base case
        if isinstance(self.value, str):
            self._format(config, skip)

    def map(self, df: pd.DataFrame) -> pd.Series:
        """Return a pandas Series by filling substitutions in `self.value`
        with values from `df`.
        """
        subs = self.find_subsitutions()
        # TODO: we should check against an "array flag
        if subs is None or any(s.root != "array_row" for s in subs):
            raise ValueError("No array substitutions to be made. Use <array_row.col_name> "
                             "to access column 'col_name' from a row of the array job DataFrame.")

        subs = set(subs)
        func_dict = {}
        template = self.value
        for i, sub in enumerate(subs):
            func_name = f"f{i}"
            template = sub.apply(template, style=func_name)
            func_dict[func_name] = parse_getter_string(sub.key)

        def apply_func(x):
            return template.format(**{k: str(v(x)) for k, v in func_dict.items()})

        return df.apply(apply_func, axis=1)


def parse_getter_string(s: str) -> Callable:
    """Take string of Python code for getting attributes and items and return
    function that does the same sequence of "getters" applied to first object.

    For instance 'obj.attr1[key1].attr2.attr3[key2]' --> function f so that
    f(obj) = obj.attr1[key1].attr2.attr3[key2].
    """
    # look for "." and "[" to know when to getattr and when to __getitem__
    # and also look for strings that are valid Python variable names
    variable_str = r"[a-zA-Z0-9_]+"
    token_specification = [
        ("getattr", r"\."),
        ("variable", variable_str),
        ("index", r"\[")
    ]
    token_pat = re.compile("|".join(f"(?P<{name}>{pat})" for name, pat in token_specification))
    tokens = list(token_pat.finditer(s))

    if not tokens[0].lastgroup == "variable":
        raise ValueError(f"Pattern {s} does not start with a valid variable name.")

    funcs = []
    # get a sequence of functions that do f(x) = x.attr or f(x) = x[key]
    #
    # we ignore the first token because it needs to be put through the function
    # that we will return
    for op, var in zip(tokens[1::2], tokens[2::2]):
        if not op.lastgroup in ["getattr", "index"] or not var.lastgroup == "variable":
            raise ValueError("Invalid pattern.")
        if op.lastgroup == "getattr":
            funcs.append(attrgetter(var[0]))
        else:
            funcs.append(itemgetter(var[0]))

    # compose all of the functions, starting with the first
    # note that the result of the function passed to the first argument
    # of reduce must be a function for this to make sense
    return reduce(lambda f, g: (lambda x: g(f(x))), funcs)


PT = TypeVar("PT", bound="Params")  # for classmethod typing


@dataclass
class Params:
    params: dict[str, Param]

    @classmethod
    def from_dict(cls: type[PT], d: dict) -> PT:
        return cls({k: Param(v, key=k) for k, v in d.items()})

    def to_dict(self) -> dict:
        return {k: v.value for k, v in self.params.items()}

    def __getitem__(self, key: str) -> Param:
        return self.params[key]

    def __iter__(self) -> Iterator:
        yield from self.params.keys()

    def format(self, config: Config, skip: list[str] | None = None) -> None:
        for v in self.params.values():
            v.format(config, skip)

    def update(self, other: Params | dict, priority: Literal["self", "other"] = "self") -> None:
        if isinstance(other, dict):
            other = Params.from_dict(other)

        if priority == "self":
            self.params = other.params | self.params
        elif priority == "other":
            self.params = self.params | other.params
        else:
            raise ValueError(f"priority must be 'self' or 'other', not '{priority}'.")


CT = TypeVar("CT", bound="Combos")  # for classmethod typing


@dataclass
class Combos:
    param_lists: dict[str, list[Param]]
    names: dict[str, list[str]] = field(default_factory=dict)
    name_template: str | None = None

    @classmethod
    def from_conf(cls: type[CT], combos_conf: dict) -> CT:
        combos_conf = combos_conf.copy()  # avoid mutating input ...probably not important here
        names = combos_conf.pop("_names", None)
        name_template = combos_conf.pop("_name_template", None)

        if isinstance(name_template, str):
            name_template = name_template.replace("^", "{").replace("$", "}")

        param_lists = {}
        for k, v in combos_conf.items():
            v = make_iterable(v)
            param_lists[k] = [Param(value=x, key=k) for x in v]

        return cls(param_lists, names, name_template)

    def __getitem__(self, key) -> dict:
        if key != "_names":
            raise KeyError(f"Can only access `_names` by key from Combos; received `{key}`.")
        return self.names

    def get_params(self) -> list[Params]:
        return [Params(x) for x in flatten(self.param_lists)]

    def get_names(self) -> list[str]:
        names_dicts = self.parse_names()
        return [self.make_name(nd) for nd in flatten(names_dicts)]

    def get_experiments(self, setup: dict, slurm: dict, dates: Dates | None) -> list[Experiment]:
        names_flat = self.get_names()
        params_flat = self.get_params()

        # fill any names from the combo names
        for params, name_dict in zip(params_flat, flatten(self.parse_names())):
            for param in params.params.values():
                if isinstance(param.value, str):
                    param.value = param.value.replace("^", "{").replace("$", "}").format(**name_dict)

        return [
            Experiment(name, params, dates=dates, setup=setup, slurm=slurm)
            for name, params in zip(names_flat, params_flat)
        ]

    # name parsing for combos
    def parse_names(self) -> dict:
        names_dict = {}
        for k, v in self.param_lists.items():
            if k in self.names:
                names_dict[k] = self.names[k]
            elif not isinstance(v, list):
                names_dict[k] = None
            else:
                names_dict[k] = list(map(str, v))
        return names_dict

    def make_name(self, names_dict: dict) -> str:
        if self.name_template is not None:
            return self.name_template.format(**names_dict)
        name_strings = []
        for k, v in names_dict.items():
            if v is None:
                continue
            if k in self.names:
                # don't add key for custom names
                name_strings.append(v)
            else:
                name_strings.append(f"{k}_{v}")

        return "_".join(name_strings)


@dataclass
class Experiment:
    """A single experiment.

    Experiment objects should be created by Config.

    They contain all of the info necessary to create a new directory, ini file,
    and slurm script for an experiment.

    Note: "experiments" in a toml config file will ultimately result in an Experiment object,
    but the Config class needs to add the general info to the experiment-specific config.
    """

    name: str
    params: Params
    dates: Optional[Dates] = None
    setup: dict = field(default_factory=dict)
    slurm: dict = field(default_factory=dict)

    def format(self, config: Config) -> None:
        self.params.format(config, skip=skip_flags)

        # format name
        tmp = Param(self.name)
        tmp.format(config)
        self.name = tmp.value

        # format setup
        tmp_setup = {}
        for k, v in self.setup.items():
            tmp = Param(v)
            tmp.format(config)
            tmp_setup[k] = tmp.value

        self.setup = tmp_setup

    def _find_array_params(self) -> list[str]:
        """Find params that have substitutions matching the array_row."""
        result = []
        for param in self.params:
            subs = self.params[param].find_subsitutions()
            if subs is not None and any(s.root == "array_row" for s in subs):
                result.append(param)
        return result

    def make_array_df(self) -> pd.DataFrame:
        array_params = self._find_array_params()
        if self.dates is None:
            if array_params:
                raise ValueError("Parameters for array job dataframe found, but no dates info provided.")
            else:
                return pd.DataFrame()

        df = self.dates.to_df()

        # make dictionary mapping param names to columns, to append to df
        arrays_dict = {}
        for name in array_params:
            param = self.params[name]
            col = param.map(df)
            arrays_dict[name] = col

        if arrays_dict:
            kwargs = (pd.DataFrame.from_dict(arrays_dict)
                   .apply(lambda x: x.to_dict(), axis=1)  # combine columns into dict because we need to pass them via --kwargs
                   .rename("kwargs")
                   )
            return df.join(kwargs)

        return df

    def make(self) -> None:
        """Make experiment directory with associated files."""
        kwargs = self.params.to_dict()
        dates_df = self.make_array_df()

        for name in self._find_array_params():
            del kwargs[name]

        ini_template_path = Path(self.setup["ini_file"])
        job_root_path = Path(self.setup["job_root"])

        python_venv = self.setup.get("python_venv", None)
        conda_venv = self.setup.get("conda_venv", None)

        # set up directory to hold ini. slurm scipt, config, results, logs
        job_name = self.setup["job_name"]

        if "out_prefix" in self.setup:
            out_name = self.setup["out_prefix"] + "_" + self.name
        else:
            out_name = self.name

        out_path = job_root_path / out_name

        if out_path.exists():
            j = 0
            while out_path.exists():
                new_out_name = out_name + str(j)
                out_path = job_root_path / new_out_name
                j += 1

        out_path.mkdir(parents=True)

        # get optional readme text
        readme = kwargs.pop("_readme", None)

        # write config file with dates
        config_path = out_path / "inversion_dates.txt"
        dates_df.to_csv(config_path, sep="\t")

        # write ini file with updated kwargs
        ini_out_path = out_path / f"{job_name}.ini"

        kwargs["outputpath"] = str(out_path)

        if "outputname" not in kwargs:
            kwargs["outputname"] = job_name

        if "merged_data_dir" not in self.setup and "merged_data_dir" not in kwargs:
            kwargs["merged_data_dir"] = str(job_root_path / "merged_data")

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
            n_array_jobs=self.dates.n_periods if self.dates else 1,
            n_kwargs=len(dates_df.columns),
            **self.slurm,
        )

        with open(out_path / "slurm.sh", "w") as f:
            f.write(slurm_script)

        with open(out_path / "readme.txt", "w") as f:
            f.write(f"Experiment: {self.name}")
            f.write("\n\n")

            if readme is not None:
                f.write(readme)
                f.write("\n\n")

            f.write("Experiment parameters:\n\n")
            for k, v in kwargs.items():
                f.write(f"{k}: {v}\n")

        print(f"Array config file and inversion wrapper for experiment {self.name} written to {out_path}.")

   
ConfT = TypeVar("ConfT", bound="Config")  # for classmethod typing


@dataclass
class Config:
    general: Params
    dates: Optional[Dates] = None
    setup: dict = field(default_factory=dict)
    slurm: dict = field(default_factory=dict)
    combos: Optional[Combos] = None
    experiments: list[Experiment] = field(default_factory=list)

    @classmethod
    def from_conf(cls: type[ConfT], toml_path: Union[str, Path], general: Optional[Params] = None) -> ConfT:
        conf = load_conf(toml_path)
        dates = Dates(**conf["dates"])
        setup = conf["setup"]
        slurm = conf["slurm"]

        if general is None:
            general = Params.from_dict(conf["general"])

        if "combos" not in conf:
            combos = None
        else:
            combos = Combos.from_conf(conf["combos"])

        if "experiments" not in conf:
            experiments = None
        else:
            experiments = [
                Experiment(name=k, params=Params.from_dict(v), dates=dates, setup=setup, slurm=slurm)
                for k, v in conf["experiments"].items()
            ]

        if experiments:
            result = cls(
                general=general, dates=dates, setup=setup, slurm=slurm, combos=combos, experiments=experiments
            )
        else:
            result = cls(general=general, dates=dates, setup=setup, slurm=slurm, combos=combos)

        # add any extra sections
        for k, v in conf.items():
            if k not in result.__dict__:
                result.__dict__[k] = v

        return result

    def __post_init__(self) -> None:
        if self.combos is not None:
            self.experiments.extend(
                self.combos.get_experiments(dates=self.dates, setup=self.setup, slurm=self.slurm)
            )
        for exp in self.experiments:
            exp.params.update(self.general)

    def format(self) -> None:
        for exp in self.experiments:
            exp.format(self)


def get_configs(toml_path: Union[str, Path]) -> list[Config]:
    """Parse config when there are multiple values for parameters in "general", specified
    via "general.combos".
    """
    conf = load_conf(toml_path)

    if "combos" in conf["general"]:
        general_combos = Combos.from_conf(conf["general"]["combos"])
        general_params = general_combos.get_params()

        conf_copy = conf["general"].copy()
        del conf_copy["combos"]
        for x in general_params:
            x.update(conf_copy)
    else:
        general_params = [conf["general"]]

    return [Config.from_conf(toml_path=toml_path, general=general) for general in general_params]
