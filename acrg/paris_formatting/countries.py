"""
Module with code related to country maps.
"""
from __future__ import annotations
from typing import cast, Optional, TypeVar, Union

from openghg_inversions import convert
import xarray as xr
from xarray.core.common import DataWithCoords

from helpers import get_area_grid, get_xr_dummies, sparse_xr_dot
from process_rhime_output import InversionOutput


# type for xr.Dataset *or* xr.DataArray
xrData = TypeVar("xrData", bound=DataWithCoords)


class Countries:
    def __init__(self, countries: xr.Dataset, country_selections: Optional[list[str]] = None) -> None:
        """Create Countries object given country map Dataset and optional list of countries to select.

        Args:
            countries: country map Dataset with `country` and `name` data variables.
            country_selections: optional list of country names to select.
        """
        self.matrix = get_xr_dummies(countries.country, cat_dim="country", categories=countries.name.values)
        self.area_grid = get_area_grid(countries.lat, countries.lon)

        if country_selections is not None:
            # check that selected countries are in the `name` variable of `countries` Dataset
            selections_check = []
            all_countries = list(map(lambda x: str(x).lower(), countries.name.values))
            for selection in country_selections:
                if selection.lower() not in all_countries:
                    selections_check.append(selection)
            if selections_check:
                raise ValueError(
                    f"Selected country/countries are not in `name` variable of `countries` Dataset: {selections_check}"
                )

            # only keep selected countries in country matrix
            filt = self.matrix.country.isin(country_selections)
            self.matrix = self.matrix.where(filt, drop=True)
            self.country_selections = country_selections
        else:
            self.country_selections = list(self.matrix.country.values)


    def get_x_to_country_mat(
        self,
        inv_out: InversionOutput,
        sparse: bool = False,
    ) -> xr.DataArray:
        """Construct a sparse matrix mapping from x sensitivities to country totals.

        Args:
            countries: xr.Dataset from country file. Must have variables: "country" with coordinate
                dimensions ("lat", "lon"), and "name" (with no coordinate dimensions).
            hbmcmc_outs: xr.Dataset from `hbmcmc_postprocessouts`.
            flux: flux used in inversion. If a constant flux was used, then `aprioriflux` from
                `hbmcmc_outs` can be used.
            area_grid: areas of each grid cell in inversion domain.
            basis_functions: xr.DataArray with coordinate dimensions ("lat", "lon") whose values assign
                grid cells to basis function boxes.
            sparse: if True, values of returned DataArray are `sparse.COO` array.

        Returns:
            xr.DataArray with coordinate dimensions ("country", "basis_region")
        """
        # compute matrix/tensor product: country_mat.T @ (area_grid * flux * basis_mat)
        # transpose doesn't need to be taken explicitly because alignment is done by dimension name
        result = sparse_xr_dot(self.matrix, self.area_grid * inv_out.flux * inv_out.basis)

        if sparse:
            return result

        # hack since `.to_numpy()` doesn't work right with sparse arrays
        return xr.apply_ufunc(lambda x: x.todense(), result)

    @staticmethod
    def _get_country_trace(
        species: str,
        x_trace: xrData,
        x_to_country: xr.DataArray,
    ) -> xrData:
        """Calculate trace(s) for total country emissions.

        Args:
            species: name of species, e.g. "co2", "ch4", "sf6", etc.
            x_trace: xr.DataArray or xr.Dataset with coordinate dimensions ("draw", "nx").
                Note: "nx" be replaced with another name, as long as the same coordinate name
                was used in `get_x_to_country_matrix`.
            x_to_country: xr.DataArray with result from `get_x_to_country_mat`

        Returns:
            xr.DataArray with coordinate dimensions ("country", "draw")

        TODO: there is a "country unit" conversion in the old code, but it seems to always product
              1.0, based on how it is used in hbmcmc
        """
        raw_trace = sparse_xr_dot(x_to_country, x_trace)
        molar_mass = convert.molar_mass(species)
        return raw_trace * 365 * 24 * 3600 * molar_mass  # type: ignore

    def get_country_trace(self, species: str, inv_out: InversionOutput) -> xr.Dataset:
        """Calculate trace(s) for total country emissions.

        Args:
            species: name of species, e.g. "co2", "ch4", "sf6", etc.
            inv_out: InversionOutput

        Returns:
            xr.Dataset with coordinate dimensions ("country", "draw")

        TODO: there is a "country unit" conversion in the old code, but it seems to always product
              1.0, based on how it is used in hbmcmc
        """
        x_to_country_mat = self.get_x_to_country_mat(inv_out)
        x_trace = inv_out.get_trace_dataset(convert_nmeasure=False, var_names="x")

        country_traces = Countries._get_country_trace(species, x_trace, x_to_country_mat)

        rename_dict = {dv: "country_" + str(dv).split("_")[1] for dv in country_traces.data_vars}
        return country_traces.rename_vars(rename_dict)

    def merge(self, other: Union[Countries, list[Countries]]) -> None:
        """Merge in another Countries object (in-place).

        Args:
            other: (list of) Countries object(s) to merge into this Country object.

        Returns:
            None (updates in-place)
        """
        if isinstance(other, Countries):
            other = [other]

        all_country_selections = self.country_selections.copy()
        for countries in other:
            all_country_selections.extend(countries.country_selections)

        if len(set(all_country_selections)) < len(all_country_selections):
            raise ValueError(
                "Duplicate countries selected. Make sure `country_selections` are disjoint."
            )

        self.country_selections = all_country_selections
        self.matrix = xr.concat([self.matrix] + [countries.matrix for countries in other], dim="country")
