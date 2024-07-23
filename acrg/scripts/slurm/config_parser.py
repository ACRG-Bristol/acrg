from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Any, Literal, Optional, TypeVar, Union

from helpers import flatten

try:
    import tomllib
except ImportError:
    import pip._vendor.tomli as tomllib


def load_conf(toml_path: Union[str, Path]) -> dict:
    with open(toml_path, "rb") as f:
        conf = tomllib.load(f)
    return conf


@dataclass
class Dates:
    year: int
    n_periods: int
    frequency: Literal["monthly", "annual"] = "annual"
    initial_month: int = 1
    array_job_id: bool = True


# regexes for toml keys and substitution parameters in config files
# The idea is that "<general.species>" will be replaced by the value of "species" in
# the "general" section of the config file.
#
# An exception: "<dates.year>" is treated specially, since array jobs may span many years
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
    def root(self) -> str:
        return self.keys[0]

    def apply(self, value: str) -> str:
        left = value[: self.start]
        right = value[self.end :]
        center = "{" + f"config.{self.root}" + "".join([f"[{key}]" for key in self.keys[1:]]) + "}"
        return left + center + right


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

    def _format(self, config: Config) -> None:
        if isinstance(self.value, str):
            # iteratively replace <key1.key2> with {config.key1[key2]} to set up for formatting below
            while s := self._find_next_sub(skip=["dates"]):
                self.value = s.apply(self.value)

            self.value = self.value.format(config=config)

    def format(self, config: Config) -> None:
        """Recursively apply _format"""
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
            self._format(config)


PT = TypeVar("PT", bound="Params")  # for classmethod typing


@dataclass
class Params:
    params: dict[str, Param]

    @classmethod
    def from_dict(cls: type[PT], d: dict) -> PT:
        return cls({k: Param(v, key=k) for k, v in d.items()})

    def __getitem__(self, key: str) -> Param:
        return self.params[key]

    def format(self, config: Config) -> None:
        for v in self.params.values():
            v.format(config)

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

    @classmethod
    def from_conf(cls: type[CT], combos_conf: dict) -> CT:
        combos_conf = combos_conf.copy()  # avoid mutating input ...probably not important here
        names = combos_conf.pop("_names", None)

        param_lists = {}
        for k, v in combos_conf.items():
            param_lists[k] = [Param(value=x, key=k) for x in v]

        return cls(param_lists, names)

    def get_params(self) -> list[Params]:
        return [Params(x) for x in flatten(self.param_lists)]

    def get_names(self) -> list[str]:
        names_dicts = self.parse_names()
        return [self.make_name(nd) for nd in flatten(names_dicts)]

    def get_experiments(self, setup: dict, slurm: dict, dates: Dates | None) -> list[Experiment]:
        names_flat = self.get_names()
        params_flat = self.get_params()

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
        self.params.format(config)

        # format name
        tmp = Param(self.name)
        tmp.format(config)
        self.name = tmp.value


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

        for x in general_params:
            x.update(conf["general"])
    else:
        general_params = [conf["general"]]

    return [Config.from_conf(toml_path=toml_path, general=general) for general in general_params]
