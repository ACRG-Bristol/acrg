#!/usr/bin/env python


get_ipython().run_line_magic("load_ext", "autoreload")
get_ipython().run_line_magic("autoreload", "2")


import numpy as np
import pandas as pd


import xarray as xr


from file_processing import get_netcdf_files, get_rhime_outs


from country_totals import get_country_trace, get_x_to_country_mat


from country_totals import get_xr_dummies, sparse_xr_dot, make_quantiles


#!pip install sparse


species = "sf6"


files = get_netcdf_files("/home/brendan/Documents/inversions/plotting/sf6_best")


outs = get_rhime_outs(files)


countries = xr.open_dataset(
    "/home/brendan/Documents/inversions/openghg_inversions/countries/country_EUROPE.nc"
)


countries_ukmo = xr.open_dataset(
    "/home/brendan/Documents/inversions/openghg_inversions/countries/country-ukmo_EUROPE.nc"
)


x_to_country_mats = [get_x_to_country_mat(countries, ds, sparse=True) for ds in outs]
x_to_country_mats_ukmo = [get_x_to_country_mat(countries_ukmo, ds, sparse=True) for ds in outs]


country_traces = [
    get_country_trace(countries, species, hbmcmc_outs=ds, x_to_country=mat)
    for ds, mat in zip(outs, x_to_country_mats)
]


country_traces_ukmo = [
    get_country_trace(countries_ukmo, species, hbmcmc_outs=ds, x_to_country=mat)
    for ds, mat in zip(outs, x_to_country_mats_ukmo)
]


country_traces = [
    trace.expand_dims({"time": [out_ds.Ytime.min().values]}) for trace, out_ds in zip(country_traces, outs)
]
country_trace_ds = xr.concat(country_traces, dim="time")


country_traces_ukmo = [
    trace.expand_dims({"time": [out_ds.Ytime.min().values]})
    for trace, out_ds in zip(country_traces_ukmo, outs)
]
country_trace_ukmo_ds = xr.concat(country_traces_ukmo, dim="time")


# # Getting country prior traces


min_model_error = 0.15


from flux_output_format import get_prior_samples


idatas = [get_prior_samples(ds, min_model_error) for ds in outs]


prior_x_traces = [idata.prior.x.isel(chain=0).rename({"x_dim_0": "basis_region"}) for idata in idatas]
country_prior_traces = [
    get_country_trace(countries, species, hbmcmc_outs=ds, x_to_country=mat, x_trace=xtrace)
    for ds, xtrace, mat in zip(outs, prior_x_traces, x_to_country_mats)
]
country_prior_traces_ukmo = [
    get_country_trace(countries_ukmo, species, hbmcmc_outs=ds, x_to_country=mat, x_trace=xtrace)
    for ds, xtrace, mat in zip(outs, prior_x_traces, x_to_country_mats_ukmo)
]


country_prior_traces = [
    trace.expand_dims({"time": [out_ds.Ytime.min().values]})
    for trace, out_ds in zip(country_prior_traces, outs)
]
country_prior_trace_ds = xr.concat(country_prior_traces, dim="time")


country_prior_traces_ukmo = [
    trace.expand_dims({"time": [out_ds.Ytime.min().values]})
    for trace, out_ds in zip(country_prior_traces_ukmo, outs)
]
country_prior_trace_ukmo_ds = xr.concat(country_prior_traces_ukmo, dim="time")


country_merged_ds = xr.merge(
    [
        country_prior_trace_ds.mean("draw").rename("countryapriori"),
        make_quantiles(country_prior_trace_ds, sample_dim="draw").rename("pcountryapriori"),
        country_trace_ds.mean("steps").rename("countryapost"),
        make_quantiles(country_trace_ds).rename("pcountryapost"),
    ]
)

country_ukmo_merged_ds = xr.merge(
    [
        country_prior_trace_ukmo_ds.mean("draw").rename("countryapriori"),
        make_quantiles(country_prior_trace_ukmo_ds, sample_dim="draw").rename("pcountryapriori"),
        country_trace_ukmo_ds.mean("steps").rename("countryapost"),
        make_quantiles(country_trace_ukmo_ds).rename("pcountryapost"),
    ]
)


# # Process flux
#
# We can't use the same method (producing flux traces) because they would be massing... 1000 samples for the EUROPE domain would be several gigabytes.
#
# Since the flux is constant on the basis regions, we can compute means and quantiles before mapping them to the original lat/lon domain.


basis_mats = [get_xr_dummies(ds.basisfunctions, cat_dim="basis_region") for ds in outs]


fluxes = [ds.fluxapriori for ds in outs]
x_traces = [ds.xtrace.rename({"nparam": "basis_region"}) for ds in outs]


x_means = [trace.mean("steps") for trace in x_traces]
x_quantiles = [make_quantiles(trace) for trace in x_traces]


flux_means = [sparse_xr_dot(flux * mat, mean) for flux, mean, mat in zip(fluxes, x_means, basis_mats)]


flux_quantiles = [
    sparse_xr_dot(flux * mat, quantiles) for flux, quantiles, mat in zip(fluxes, x_quantiles, basis_mats)
]


prior_x_means = [trace.mean("draw") for trace in prior_x_traces]
prior_x_quantiles = [make_quantiles(trace, sample_dim="draw") for trace in prior_x_traces]


prior_flux_means = [
    sparse_xr_dot(flux * mat, mean) for flux, mean, mat in zip(fluxes, prior_x_means, basis_mats)
]


prior_flux_quantiles = [
    sparse_xr_dot(flux * mat, quantiles)
    for flux, quantiles, mat in zip(fluxes, prior_x_quantiles, basis_mats)
]


times = [ds.Ytime.min().values for ds in outs]


prior_flux_mean_ds = xr.concat(
    [da.expand_dims({"time": [time]}) for da, time in zip(prior_flux_means, times)], dim="time"
)
flux_mean_ds = xr.concat(
    [da.expand_dims({"time": [time]}) for da, time in zip(flux_means, times)], dim="time"
)


prior_flux_quantile_ds = xr.concat(
    [da.expand_dims({"time": [time]}) for da, time in zip(prior_flux_quantiles, times)], dim="time"
)
flux_quantile_ds = xr.concat(
    [da.expand_dims({"time": [time]}) for da, time in zip(flux_quantiles, times)], dim="time"
)


flux_merged_ds = xr.merge(
    [
        prior_flux_mean_ds.rename("fluxapriori"),
        prior_flux_quantile_ds.rename("pfluxapriori"),
        flux_mean_ds.rename("fluxapost"),
        flux_quantile_ds.rename("pfluxapot"),
    ]
)


# # Merge all


paris_countries = [
    "IRELAND",
    "UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN IRELAND",
    "FRANCE",
    "BELGIUM",
    "NETHERLANDS",
    "GERMANY",
    "DENMARK",
    "SWITZERLAND",
    "AUSTRIA",
    "ITALY",
    "CZECHIA",
    "POLAND",
    "HUNGARY",
    "SLOVAKIA",
    "NORWAY",
    "SWEDEN",
    "FINLAND",
]


country_filt = country_merged_ds.country.isin(paris_countries)


country_ukmo_filt = country_ukmo_merged_ds.country.isin(["BENELUX", "RestEU", "SpaPor"])


country_merged_ds.where(country_filt, drop=True)


country_final_ds = xr.concat(
    [
        country_merged_ds.where(country_filt, drop=True),
        country_ukmo_merged_ds.where(country_ukmo_filt, drop=True),
    ],
    dim="ncountries",
)


rhime_emissions = xr.merge([flux_merged_ds, country_final_ds]).squeeze()


rhime_emissions
