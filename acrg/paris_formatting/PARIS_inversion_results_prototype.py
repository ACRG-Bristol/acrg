#!/usr/bin/env python
# coding: utf-8

# ## PARIS inverse modelling results
#
# This notebook contains functions that can be used to plot and compare results from different inverse model outputs. This includes options to compare results between models and compare results to apriori estimates:
#
# - Posterior country fluxes, total from all sectors (for ELRIS, InTEM and RHIME)
# - Posterior modelled total mole fractions (currently just InTEM - will be updated to include other models soon)
# - Posterior modelled baseline mole fractions (currently just InTEM - will be updated to include other models soon)
#
# Future updates may include:
#
# - Comparisons between model and inventory emissions estimates
# - Sector-level emissions
# - Comparison between each model's country/region definition

# ### Notebook setup:
#
# Run the three cells below, before running any of the plotting code.
#
# Edit the data path to point towards where the model output is.


data_dir = "/group/chemistry/acrg/PARIS_results_sharing/"


import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import cartopy

from matplotlib.dates import YearLocator, MonthLocator
from matplotlib.ticker import NullFormatter


# Species and country options:
species_all = ["ch4", "hfc134a", "pfc218", "sf6"]
countries_all = [
    "IRELAND",
    "UK",
    "FRANCE",
    "NETHERLANDS",
    "GERMANY",
    "DENMARK",
    "SWITZERLAND",
    "AUSTRIA",
    "ITALY",
    "BELUX",
    "BENELUX",
    "NW_EU",
    "NW_EU2",
    "NW_EU_CONTINENT",
    "CW_EU",
    "EU_GRP2",
]
countrycodes_all = ["IE", "GB", "FR", "NL", "DE", "DK", "CH", "AT", "IT", "BE-LU", "BE-LE-NE"]

countrycodes_dict = dict(zip(countries_all, countrycodes_all))

# Species and site info, used to customise plots

species_print = {
    "ch4": "CH$_4$",
    "hfc134a": "HFC-134a",
    "pfc218": "PFC-218",
    "sf6": "SF$_6$",
    "n2o": "N$_2$O",
}

intem_species = {"ch4": "ch4", "hfc134a": "h134a", "pfc218": "pfc218", "sf6": "sf6", "n2o": "n2o"}

period = {
    "intem": {"ch4": "monthly", "hfc134a": "yearly", "pfc218": "yearly", "sf6": "monthly", "n2o": "monthly"},
    "rhime": {"ch4": "monthly", "hfc134a": "yearly", "pfc218": "yearly", "sf6": "monthly", "n2o": "monthly"},
    "elris_name": {
        "ch4": "monthly",
        "hfc134a": "yearly",
        "pfc218": "yearly",
        "sf6": "yearly",
        "n2o": "monthly",
    },
    "elris_flexpart": {
        "ch4": "monthly",
        "hfc134a": "yearly",
        "pfc218": "yearly",
        "sf6": "yearly",
        "n2o": "monthly",
    },
}

units_scaling = {
    "intem": {"ch4": 1e9, "hfc134a": 1e6, "pfc218": 1e6, "sf6": 1e6, "n2o": 1e6},
    "rhime": {"ch4": 1e12, "hfc134a": 1e9, "pfc218": 1e9, "sf6": 1e9, "n2o": 1e9},
    "elris_name": {"ch4": 1e9, "hfc134a": 1e6, "pfc218": 1e6, "sf6": 1e6, "n2o": 1e6},
    "elris_flexpart": {"ch4": 1e9, "hfc134a": 1e6, "pfc218": 1e6, "sf6": 1e6, "n2o": 1e6},
}

units_print = {"ch4": "T", "hfc134a": "G", "pfc218": "G", "sf6": "G", "n2o": "G"}

dt_units = {
    "ch4": "datetime64[M]",
    "hfc134a": "datetime64[Y]",
    "pfc218": "datetime64[Y]",
    "sf6": "datetime64[M]",
    "n2o": "datetime64[M]",
}

dt_units = {
    "intem": {
        "ch4": "datetime64[M]",
        "hfc134a": "datetime64[Y]",
        "pfc218": "datetime64[Y]",
        "sf6": "datetime64[M]",
        "n2o": "datetime64[M]",
    },
    "rhime": {
        "ch4": "datetime64[M]",
        "hfc134a": "datetime64[Y]",
        "pfc218": "datetime64[Y]",
        "sf6": "datetime64[M]",
        "n2o": "datetime64[M]",
    },
    "elris_name": {
        "ch4": "datetime64[M]",
        "hfc134a": "datetime64[Y]",
        "pfc218": "datetime64[Y]",
        "sf6": "datetime64[Y]",
        "n2o": "datetime64[M]",
    },
    "elris_flexpart": {
        "ch4": "datetime64[M]",
        "hfc134a": "datetime64[Y]",
        "pfc218": "datetime64[Y]",
        "sf6": "datetime64[Y]",
        "n2o": "datetime64[M]",
    },
}


mf_units_scaling = {"ch4": 1e-9, "hfc134a": 1e-12, "pfc218": 1e-12, "sf6": 1e-12, "n2o": 1e-12}

mf_units_print = {"ch4": "ppb", "hfc134a": "ppt", "pfc218": "ppt", "sf6": "ppt", "n2o": "ppt"}

intem_site = {
    "MHD": "MH",
    "TAC": "T2",
    "CBW": "CB",
    "HFD": "H1",
    "LUT": "LW",
    "JFJ": "J1",
    "RGL": "R1",
    "BSD": "B2",
    "WAO": "WB",
    "CMN": "M5",
    "ZEP": "ZE",
    "BIR": "BI",
    "GAT": "GW",
    "HEI": "HB",
    "HEL": "HN",
    "HPB": "HP",
    "HTM": "HY",
    "HUN": "HG",
    "KIT": "KA",
    "KRE": "KR",
    "LIN": "LI",
    "NOR": "NO",
    "OPE": "OP",
    "OXK": "OK",
    "PAL": "PA",
    "SAC": "SA",
    "SSL": "SL",
    "STE": "SK",
    "TOH": "TO",
    "TRN": "TR",
    "UTO": "UT",
    "WES": "WE",
    "TTA": "AN",
}

model_labels = {
    "intem": "InTEM",
    "rhime": "RHIME",
    "elris_name": "ELRIS_NAME",
    "elris_flexpart": "ELRIS_FLEXPART",
}

model_colors = {
    "intem": ["darkslateblue", "dodgerblue"],
    "rhime": ["darkgreen", "green"],
    "elris_name": ["purple", "mediumpurple"],
    "elris_flexpart": ["darkorange", "goldenrod"],
}

model_fp = {
    "intem": "InTEM_NAME_EUROPE",
    "rhime": "RHIME_NAME_EUROPE",
    "elris_name": "ELRIS_NAME_EUROPE_EDGAR_NAME_daily",
    "elris_flexpart": "ELRIS_FLEXPART-IFS_EUROPE_EDGAR_PARIS-FP_daily",
}

model_species = {
    "intem": {"ch4": "ch4", "hfc134a": "h134a", "pfc218": "pfc218", "sf6": "sf6", "n2o": "n2o"},
    "rhime": {"ch4": "ch4", "hfc134a": "hfc134a", "pfc218": "pfc218", "sf6": "sf6", "n2o": "n2o"},
    "elris_name": {"ch4": "CH4", "hfc134a": "HFC_134a", "pfc218": "PFC_218", "sf6": "SF6", "n2o": "N2O"},
    "elris_flexpart": {"ch4": "CH4", "hfc134a": "HFC_134a", "pfc218": "PFC_218", "sf6": "SF6", "n2o": "N2O"},
}

model_q_indices = {"intem": [0, 1], "rhime": [1, 2], "elris_name": [0, 1], "elris_flexpart": [0, 1]}


font = {"size": 12}
plt.rc("font", **font)


# ### 1. Plot timeseries of country/region fluxes
#
# Edit and run this cell to choose inputs:


species = "hfc134a"
plot_regions = [
    "IRELAND",
    "UK",
    "FRANCE",
    "BENELUX",
    "GERMANY",
    "NETHERLANDS",
]  # works best with 6 countries but can run with any number (will always create at least 4 subplots though)
models = ["intem", "rhime", "elris_name"]
start_date = "2018-01-01"
end_date = "2022-01-01"
plot_inventory = False


# Then run this to plot:


ds_all = {}

for m in models:
    try:
        with xr.open_dataset(
            os.path.join(data_dir, f"{model_fp[m]}_{model_species[m][species]}_{period[m][species]}.nc")
        ) as in_ds:
            ds_all[m] = in_ds
    except:
        print(f"Cannot find {m} file for {species}.")

n_cols = int(len(plot_regions) / 2)

if n_cols <= 1:
    n_cols = 2

fig, ax = plt.subplots(2, n_cols, figsize=(n_cols * 6, 8), constrained_layout=True)

a, b = 0, 0
max_cf = []

for i, country in enumerate(plot_regions):
    if plot_inventory == True:
        try:
            with xr.open_dataset(
                os.path.join(data_dir, f'UNFCCC_inventory_{model_species["intem"][species]}.nc')
            ) as inv_ds:
                inv_c_index = np.where(inv_ds["country"].values == country)[0][0]
                ax[a, b].bar(
                    inv_ds.time.values,
                    inv_ds["inventory"].values[:, inv_c_index] / units_scaling["intem"][species],
                    np.timedelta64(340, "D"),
                    color="white",
                    edgecolor="black",
                    align="edge",
                    label="Inventory 2023",
                    zorder=0,
                )
        except:
            print(f"No inventory data available for {country}")

    for m in models:
        try:
            if m == "intem":
                c_key = "countrynames"
                country_search = country
            elif m == "rhime":
                c_key = "country"
                if country == "UK":
                    country_search = "UNITED KINGDOM"
                else:
                    country_search = country
            elif m == "elris_name" or m == "elris_flexpart":
                c_key = "country"
                country_search = countrycodes_dict[country]

            country_index = np.where(ds_all[m][c_key].values.astype(str) == country_search)[0][0]
            t0 = np.where(
                ds_all[m].time.values.astype(dt_units[m][species])
                == np.datetime64(start_date).astype(dt_units[m][species])
            )[0][0]
            t1 = (
                np.max(
                    np.where(
                        ds_all[m].time.values.astype(dt_units[m][species])
                        <= np.datetime64(end_date).astype((dt_units[m][species]))
                    )
                )
                + 1
            )
            ax[a, b].plot(
                ds_all[m].time.values[t0:t1].astype(dt_units[m][species]),
                ds_all[m]["countryapost"].values[t0:t1, country_index] / units_scaling[m][species],
                label=model_labels[m],
                color=model_colors[m][0],
            )
            ax[a, b].fill_between(
                ds_all[m].time.values[t0:t1].astype(dt_units[m][species]),
                ds_all[m]["qcountryapost"].values[t0:t1, model_q_indices[m][0], country_index]
                / units_scaling[m][species],
                ds_all[m]["qcountryapost"].values[t0:t1, model_q_indices[m][1], country_index]
                / units_scaling[m][species],
                alpha=0.3,
                color=model_colors[m][0],
            )
            ax[a, b].plot(
                ds_all[m].time.values[t0:t1].astype(dt_units[m][species]),
                ds_all[m]["countryapriori"].values[t0:t1, country_index] / units_scaling[m][species],
                label=f"{model_labels[m]} prior",
                color=model_colors[m][0],
                linestyle="dashed",
            )
            max_cf.append(
                np.nanmax(
                    (
                        np.nanmax(
                            ds_all[m]["countryapost"].values[t0:t1, country_index] / units_scaling[m][species]
                        ),
                        np.nanmax(
                            ds_all[m]["countryapriori"].values[t0:t1, country_index]
                            / units_scaling[m][species]
                        ),
                    )
                )
            )
        except:
            print(
                f"ERROR: Either start and end dates are incorrect or there is no {country} emissions in {m}."
            )
            print(f"Skipping plotting {m}.")

    # format each subplot
    ax[a, b].set_ylabel(f"{species_print[species]} ({units_print[species]}g y$^{{-1}}$)")
    ax[a, b].set_xlim(
        [
            np.datetime64(start_date).astype("datetime64[M]") - np.timedelta64(1, "M"),
            np.datetime64(end_date).astype("datetime64[M]") + np.timedelta64(1, "M"),
        ]
    )

    if len(models) == 3:
        if plot_inventory == True:
            ncol = 2
        else:
            ncol = 3
    else:
        ncol = 2

    leg = ax[a, b].legend(ncol=ncol, borderpad=0.4, columnspacing=1.0, fontsize=10)
    if plot_inventory == True:
        for l in leg.legend_handles[:-1]:
            l.set_linewidth(3.0)
    else:
        for l in leg.legend_handles:
            l.set_linewidth(3.0)

    ax[a, b].set_title(f"{country}")
    ax[a, b].grid(visible=True, which="major", alpha=0.4)
    ax[a, b].xaxis.set_minor_locator(MonthLocator())
    ax[a, b].xaxis.set_minor_formatter(NullFormatter())
    ax[a, b].xaxis.set_major_locator(YearLocator())

    # increase row and column counts
    if (b - (n_cols - 1)) == 0:
        b = 0
        a += 1
    else:
        b += 1

for a in range(2):
    for b in range(n_cols):
        # ax[a,b].set_ylim([0,(np.max(max_cf)+(0.1*np.max(max_cf)))])
        ax[a, b].set_ylim(bottom=0)


print("NOTE: If all the data is not within axis limits, adjust the set_ylim parameter")


# Save plot here:


# output_path = '/home/h02/aramsden/results/PARIS_results_comparison/n2o_countryflux_comparison_IE_UK_no_ylim'

# fig.savefig(output_path,bbox_inches='tight',pad_inches=0.2,dpi=300)


# ### 2. Plot timeseries of modelled and observed mole fractions
#
# Edit and run this cell to choose inputs:


species = "hfc134a"
site = "MHD"
models = ["intem", "elris_name"]
start_date = "2018-01-01"
end_date = "2023-01-01"


# Then run this cell to plot:


ds_all = {}

for m in models:
    try:
        ds_all[m] = xr.open_dataset(
            os.path.join(
                data_dir, f"{model_fp[m]}_{model_species[m][species]}_{period[m][species]}_concentrations.nc"
            )
        )
    except:
        print(f"Cannot find {m} file for {species}.")

min_mf = []
max_mf = []
ax_all = []
ax2_all = []
annotate_coords1 = (0.1, 0.8)
annotate_coords2 = (0.1, 0.6)

fig = plt.figure(constrained_layout=True, figsize=(15, len(models) * 3))
gs = fig.add_gridspec(len(models), 2, width_ratios=[0.8, 0.2])

for i, m in enumerate(models):
    if m == "intem":
        search_site = intem_site[site]
    else:
        search_site = site

    # try:

    t0 = np.min(np.where(ds_all[m].time.values >= np.datetime64(start_date).astype("datetime64[ns]")))
    t1 = np.max(np.where(ds_all[m].time.values <= np.datetime64(end_date).astype("datetime64[ns]"))) + 1
    s = np.where(ds_all[m]["sitenames"].values.astype("str") == search_site)[0][0]

    ax = fig.add_subplot(gs[i, 0])
    ax2 = fig.add_subplot(gs[i, 1])
    ax_all.append(ax)
    ax2_all.append(ax2)

    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["Yapriori"].values[t0:t1, s] / mf_units_scaling[species],
        color=model_colors[m][1],
        label=f"{model_labels[m]} apriori mf",
    )
    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["Yapost"].values[t0:t1, s] / mf_units_scaling[species],
        color=model_colors[m][0],
        label=f"{model_labels[m]} posterior mean mf",
    )
    ax.fill_between(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["qYapost"].values[t0:t1, 0, s] / mf_units_scaling[species],
        ds_all[m]["qYapost"].values[t0:t1, 1, s] / mf_units_scaling[species],
        color=model_colors[m][0],
        alpha=0.2,
        label=f"{model_labels[m]} posterior std dev mf",
    )
    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species],
        color="black",
        label=f"Obs ({model_labels[m]})",
    )

    min_mf.append(np.nanmin(ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]))
    max_mf.append(np.nanmax(ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]))

    ax.set_title("InTEM")
    ax.set_ylabel(f"{species_print[species]} {site} ({mf_units_print[species]})")
    leg = ax.legend(ncol=2, borderpad=0.2, columnspacing=1.0)
    for l in leg.legend_handles:
        l.set_linewidth(3.0)

    if (
        int(
            np.datetime64(end_date).astype("datetime64[M]")
            - np.datetime64(start_date).astype("datetime64[M]")
        )
        > 12
    ):
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(YearLocator())
    else:
        ax.xaxis.set_major_locator(MonthLocator())

    diff = (
        ds_all[m]["Yapost"].values[t0:t1, s] / mf_units_scaling[species]
        - ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]
    )
    diff_mean = np.round(np.nanmean(diff), 2)
    diff_sd = np.round(np.nanstd(diff), 2)

    pr_diff = (
        ds_all[m]["Yapriori"].values[t0:t1, s] / mf_units_scaling[species]
        - ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]
    )
    pr_diff_mean = np.round(np.nanmean(pr_diff), 2)
    pr_diff_sd = np.round(np.nanstd(pr_diff), 2)

    a, b, c = ax2.hist(diff, bins=30, color=model_colors[m][0])
    a_pr, b_pr, c_pr = ax2.hist(pr_diff, bins=30, color=model_colors[m][1], alpha=0.7)

    with np.printoptions(precision=2, suppress=True):
        ax2.annotate(
            "$\mu$: " + str(diff_mean) + "\n$\sigma$: " + str(diff_sd),
            xy=annotate_coords1,
            xycoords="axes fraction",
            color=model_colors[m][0],
        )
        ax2.annotate(
            "$\mu$: " + str(pr_diff_mean) + "\n$\sigma$: " + str(pr_diff_sd),
            xy=annotate_coords2,
            xycoords="axes fraction",
            color=model_colors[m][1],
        )
        ax2.vlines(0, 0, np.max(a), color="black", linewidth=3.0)

    # except:
    #    print(f'ERROR: No site {site} in model {m} for this time period.')
    #    print(f'Skipping plotting {m}.')


for i in range(len(models)):
    ax_all[i].set_ylim([min(min_mf) - (0.02 * min(min_mf)), max(max_mf) + (0.05 * max(max_mf))])

print("NOTE: If all the data is not within axis limits, adjust the set_ylim")
print("NOTE: If annotations in the histograms are not displaying correctly, adjust annotate_coords.")


# Save plot:


output_path = None

# fig.savefig(output_path,bbox_inches='tight',pad_inches=0.2,dpi=300)


# ### 3. Plot timeseries of modelled baselines
#
# This uses the same input settings as the mole fraction plot.
#
# Run this cell to plot:


ds_all = {}

for m in models:
    try:
        ds_all[m] = xr.open_dataset(
            os.path.join(
                data_dir, f"{model_fp[m]}_{model_species[m][species]}_{period[m][species]}_concentrations.nc"
            )
        )
    except:
        print(f"Cannot find {m} file for {species}.")

min_mf = []
max_mf = []
ax_all = []
ax2_all = []
annotate_coords1 = (0.1, 0.8)
annotate_coords2 = (0.1, 0.6)

fig = plt.figure(constrained_layout=True, figsize=(15, len(models) * 3))
gs = fig.add_gridspec(len(models), 2, width_ratios=[0.8, 0.2])

for i, m in enumerate(models):
    if m == "intem":
        search_site = intem_site[site]
    else:
        search_site = site

    # try:

    t0 = np.min(np.where(ds_all[m].time.values >= np.datetime64(start_date).astype("datetime64[ns]")))
    t1 = np.max(np.where(ds_all[m].time.values <= np.datetime64(end_date).astype("datetime64[ns]"))) + 1
    s = np.where(ds_all[m]["sitenames"].values.astype(str) == search_site)[0][0]

    ax = fig.add_subplot(gs[i, 0])
    ax2 = fig.add_subplot(gs[i, 1])
    ax_all.append(ax)
    ax2_all.append(ax2)

    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species],
        color="grey",
        label=f"Obs ({model_labels[m]})",
    )
    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["YaprioriBC"].values[t0:t1, s] / mf_units_scaling[species],
        color=model_colors[m][1],
        label=f"{model_labels[m]} apriori baseline",
    )
    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["YapostBC"].values[t0:t1, s] / mf_units_scaling[species],
        color=model_colors[m][0],
        label=f"{model_labels[m]} posterior mean baseline",
    )

    # if m != 'intem':
    #    ax.fill_between(ds_all[m].time.values[t0:t1],
    #                    ds_all[m]['qYapostBC'].values[t0:t1,0,s],
    #                    ds_all[m]['qYapostBC'].values[t0:t1,1,s],
    #                    color=model_colors[m][0],alpha=0.2,label=f'{model_labels[m]} posterior std dev baseline')

    min_mf.append(np.nanmin(ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]))
    max_mf.append(np.nanmax(ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]))

    ax.set_title("InTEM")
    ax.set_ylabel(f"{species_print[species]} {site} ({mf_units_print[species]})")
    leg = ax.legend(ncol=4, borderpad=0.2, columnspacing=1.0)
    for l in leg.legend_handles:
        l.set_linewidth(3.0)

    if (
        int(
            np.datetime64(end_date).astype("datetime64[M]")
            - np.datetime64(start_date).astype("datetime64[M]")
        )
        > 12
    ):
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(YearLocator())
    else:
        ax.xaxis.set_major_locator(MonthLocator())

    diff = (
        ds_all[m]["YapostBC"].values[t0:t1, s] / mf_units_scaling[species]
        - ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]
    )
    diff_mean = np.round(np.nanmean(diff), 2)
    diff_sd = np.round(np.nanstd(diff), 2)

    pr_diff = (
        ds_all[m]["YaprioriBC"].values[t0:t1, s] / mf_units_scaling[species]
        - ds_all[m]["Yobs"].values[t0:t1, s] / mf_units_scaling[species]
    )
    pr_diff_mean = np.round(np.nanmean(pr_diff), 2)
    pr_diff_sd = np.round(np.nanstd(pr_diff), 2)

    a, b, c = ax2.hist(diff, bins=30, color=model_colors[m][0])
    a_pr, b_pr, c_pr = ax2.hist(pr_diff, bins=30, color=model_colors[m][1], alpha=0.7)

    ax2.annotate(
        "$\mu$: " + str(diff_mean) + "\n$\sigma$: " + str(diff_sd),
        xy=annotate_coords1,
        xycoords="axes fraction",
        color=model_colors[m][0],
    )
    ax2.annotate(
        "$\mu$: " + str(pr_diff_mean) + "\n$\sigma$: " + str(pr_diff_sd),
        xy=annotate_coords2,
        xycoords="axes fraction",
        color=model_colors[m][1],
    )
    ax2.vlines(0, 0, np.max(a), color="grey", linewidth=3.0)

    # except:
    #    print(f'ERROR: No site {site} in model {m} for this time period.')
    #    print(f'Skipping plotting {m}.')

for i in range(len(models)):
    ax_all[i].set_ylim([min(min_mf) - (0.02 * min(min_mf)), max(max_mf) + (0.05 * max(max_mf))])

print("NOTE: If all the data is not within axis limits, adjust the ylim parameter")
print("NOTE: If annotations in the histograms are not displaying correctly, adjust annotate_coords.")


# Save plot:


output_path = None

# fig.savefig(output_path,bbox_inches='tight',pad_inches=0.2,dpi=300)


# ### 4. Directly compare posterior mole fractions
#
# Edit this cell to choose inputs:


species = "ch4"
site = "MHD"
models = ["intem"]  # ['intem','rhime','empa']
start_date = "2018-01-01"
end_date = "2023-01-01"


# Then run this cell to plot:


min_mf = []
max_mf = []
ax_all = []
ax2_all = []
annotate_coords1 = (0.1, 0.8)
annotate_coords2 = (0.1, 0.6)

for m in models:
    ds_all[m] = xr.open_dataset(
        os.path.join(
            data_dir, f"{model_fp[m]}_{model_species[m][species]}_{period[m][species]}_concentrations.nc"
        )
    )

fig = plt.figure(constrained_layout=True, figsize=(15, 4))
gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 0.2])

ax = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

for i, m in enumerate(models):
    # try:

    t0 = np.min(np.where(ds_all[m].time.values >= np.datetime64(start_date).astype("datetime64[ns]")))
    t1 = np.max(np.where(ds_all[m].time.values <= np.datetime64(end_date).astype("datetime64[ns]"))) + 1
    s = np.where(ds_all[m]["sitenames"].values == intem_site[site])[0][0]

    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["Yapost"].values[t0:t1, s],
        color=model_colors[m][0],
        label=f"{model_labels[m]} posterior mean mf",
    )

    ax.fill_between(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["qYapost"].values[t0:t1, 0, s],
        ds_all[m]["qYapost"].values[t0:t1, 1, s],
        color=model_colors[m][0],
        alpha=0.2,
        label=f"{model_labels[m]} posterior std dev mf",
    )

    min_mf.append(np.nanmin(ds_all[m]["Yobs"].values[t0:t1, s]))
    max_mf.append(np.nanmax(ds_all[m]["Yobs"].values[t0:t1, s]))

    ax.set_ylabel(f"{species_print[species]} {site} ({mf_units_print[species]})")

    mf_mean = np.round(np.nanmean(ds_all[m]["Yapost"].values[t0:t1, s]), 2)
    mf_sd = np.round(np.nanstd(ds_all[m]["Yapost"].values[t0:t1, s]), 2)

    a, b, c = ax2.hist(ds_all[m]["Yapost"].values[t0:t1, s], bins=30, color=model_colors[m][0])

    ax2.annotate(
        f"$\mu$: {mf_mean}\n$\sigma$: {mf_sd}",
        xy=annotate_coords1,
        xycoords="axes fraction",
        color=model_colors[m][0],
    )

    # except:
    #    print(f'ERROR: No site {site} in model {m} for this time period.')
    #    print(f'Skipping plotting {m}.')

    leg = ax.legend(ncol=2, borderpad=0.2, columnspacing=1.0)
    for l in leg.legend_handles:
        l.set_linewidth(3.0)

    if (
        int(
            np.datetime64(end_date).astype("datetime64[M]")
            - np.datetime64(start_date).astype("datetime64[M]")
        )
        > 12
    ):
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(YearLocator())
    else:
        ax.xaxis.set_major_locator(MonthLocator())

        ax.set_ylim([min(min_mf) - (0.02 * min(min_mf)), max(max_mf) + (0.05 * max(max_mf))])

print("NOTE: If all the data is not within axis limits, adjust the ylim parameter")
print("NOTE: If annotations in the histograms are not displaying correctly, adjust annotate_coords.")


# Save plot:


output_path = None

# fig.savefig(output_path,bbox_inches='tight',pad_inches=0.2,dpi=300)


# ### 5. Directly compare posterior baselines
#
# Edit this cell to choose inputs:


species = "ch4"
site = "MHD"
models = ["intem"]  # ['intem','rhime','empa']
start_date = "2018-01-01"
end_date = "2023-01-01"


# Run this cell to plot:


min_mf = []
max_mf = []
ax_all = []
ax2_all = []
annotate_coords1 = (0.1, 0.8)
annotate_coords2 = (0.1, 0.6)

for m in models:
    ds_all[m] = xr.open_dataset(
        os.path.join(
            data_dir, f"{model_fp[m]}_{model_species[m][species]}_{period[m][species]}_concentrations.nc"
        )
    )

fig = plt.figure(constrained_layout=True, figsize=(15, 4))
gs = fig.add_gridspec(1, 2, width_ratios=[0.8, 0.2])

ax = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

for i, m in enumerate(models):
    # try:

    t0 = np.min(np.where(ds_all[m].time.values >= np.datetime64(start_date).astype("datetime64[ns]")))
    t1 = np.max(np.where(ds_all[m].time.values <= np.datetime64(end_date).astype("datetime64[ns]"))) + 1
    s = np.where(ds_all[m]["sitenames"].values == intem_site[site])[0][0]

    ax.plot(
        ds_all[m].time.values[t0:t1],
        ds_all[m]["YapostBC"].values[t0:t1, s],
        color=model_colors[m][0],
        label=f"{model_labels[m]} posterior mean baseline",
    )

    if m != "intem":
        ax.fill_between(
            ds_all[m].time.values[t0:t1],
            ds_all[m]["qYapostBC"].values[t0:t1, 0, s],
            ds_all[m]["qYapostBC"].values[t0:t1, 1, s],
            color=model_colors[m][0],
            alpha=0.2,
            label=f"{model_labels[m]} posterior std dev baseline",
        )

    min_mf.append(np.nanmin(ds_all[m]["Yobs"].values[t0:t1, s]))
    max_mf.append(np.nanmax(ds_all[m]["Yobs"].values[t0:t1, s]))

    ax.set_ylabel(f"{species_print[species]} {site} ({mf_units_print[species]})")

    mf_mean = np.round(np.nanmean(ds_all[m]["YapostBC"].values[t0:t1, s]), 2)
    mf_sd = np.round(np.nanstd(ds_all[m]["Yapost"].values[t0:t1, s]), 2)

    a, b, c = ax2.hist(ds_all[m]["YapostBC"].values[t0:t1, s], bins=30, color=model_colors[m][0])

    ax2.annotate(
        f"$\mu$: {mf_mean}\n$\sigma$: {mf_sd}",
        xy=annotate_coords1,
        xycoords="axes fraction",
        color=model_colors[m][0],
    )

    # except:
    #    print(f'ERROR: No site {site} in model {m} for this time period.')
    #    print(f'Skipping plotting {m}.')

    leg = ax.legend(ncol=2, borderpad=0.2, columnspacing=1.0)
    for l in leg.legend_handles:
        l.set_linewidth(3.0)

    if (
        int(
            np.datetime64(end_date).astype("datetime64[M]")
            - np.datetime64(start_date).astype("datetime64[M]")
        )
        > 12
    ):
        ax.xaxis.set_minor_locator(MonthLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_major_locator(YearLocator())
    else:
        ax.xaxis.set_major_locator(MonthLocator())

        ax.set_ylim([min(min_mf) - (0.02 * min(min_mf)), max(max_mf) + (0.05 * max(max_mf))])

print("NOTE: If all the data is not within axis limits, adjust the ylim parameter")
print("NOTE: If annotations in the histograms are not displaying correctly, adjust annotate_coords.")


# Save plot:


output_path = None

# fig.savefig(output_path,bbox_inches='tight',pad_inches=0.2,dpi=300)


# ### 6. Plot posterior country fluxes - lat lon grid comparison
#
# THIS NEEDS UPDATING TO DIRECTLY COMPARE MULTIPLE MODELS - CURRENTLY THESE NEED TO BE INPUT MANUALLY


species = "hfc134a"  # select the species you want to plot
date = "2021-01-01"  # select the time period to plot (a month or year, depending on the inversion period)

plot_type = "posterior"  # posterior or prior
plot_area = "CWEU"  # options for: UK, FRANCE, GERMANY, NWEU, CWEU


# Read in emissions output

intem = xr.open_dataset(
    os.path.join(data_dir, f"InTEM_NAME_EUROPE_{intem_species[species]}_{period[m][species]}.nc")
)
# rhime = xr.open_dataset(os.path.join(data_dir,f'RHIME_NAME_EUROPE_{species}_{period}.nc'))
# empa = xr.open_dataset(os.path.join(data_dir,f'{empa_name}_FLEXPART_EUROPE_{species}_{period}.nc'))

plot_lat = np.hstack((intem.lat.values, intem.lat.values[-1] + (intem.lat.values[1] - intem.lon.values[0])))

try:
    intem_t0 = np.where(intem.time.values == np.datetime64(date).astype("datetime64[ns]"))[0][0]
except:
    print(f"No InTEM emissions estimate for {date}")

fluxlim = {"ch4": [0, 1e-7], "hfc134a": [0, 1e-11], "pfc218": [0, 1e-14], "sf6": [0, 5e-13], "n20": [0, 1e-8]}

difflim = {
    "ch4": [-1e-7, 1e-7],
    "hfc134a": [-1e-11, 1e-11],
    "pfc218": [-1e-14, 1e-14],
    "sf6": [-5e-13, 5e-13],
    "n20": [-1e-8, 1e-8],
}

region_limits = {
    "UK": [-12, 4, 49, 62],  # min_lon, max_lon, min_lat, max_lat
    "FRANCE": [-6, 9, 42, 52],
    "GERMANY": [2, 18, 45, 60],
    "NWEU": [-11, 11, 45, 62],
    "CWEU": [-12, 27, 37, 66],
}


cmap = "viridis"
cmap_diff = "bwr"

if plot_type == "posterior":
    ds_key = "fluxapost"
    print_type = f"Posterior mean"
elif plot_type == "prior":
    ds_key = "fluxapriori"
    print_type = f"Prior mean"

fig, ax = plt.subplots(
    2, 3, figsize=(15, 6), constrained_layout=True, subplot_kw={"projection": cartopy.crs.PlateCarree()}
)

for i in range(2):
    for j in range(3):
        ax[i, j].add_feature(cartopy.feature.BORDERS, linestyle=":", edgecolor="beige", linewidth=1.0)
        ax[i, j].coastlines(resolution="50m", color="beige", linewidth=1.0)
        ax[i, j].set_extent(region_limits[plot_area])

# plot intem
ax[0, 0].pcolormesh(
    intem.lon.values,
    intem.lat.values,
    intem[ds_key][intem_t0, :-1, :-1],
    cmap="viridis",
    vmin=fluxlim[species][0],
    vmax=fluxlim[species][1],
    shading="flat",
)

ax[0, 0].set_title("InTEM")

# flux colorbar
levels = np.linspace(fluxlim[species][0], fluxlim[species][1])
m = plt.cm.ScalarMappable(cmap=cmap)
m.set_array(levels)
m.set_clim(fluxlim[species])

color_bar1 = fig.colorbar(
    m, orientation="vertical", cmap=cmap, extend="max", ax=ax[0, :], shrink=0.9, pad=0.005
)
color_bar1.set_label(f"{print_type}\n{species_print[species]} flux {date}\n(mol m$^{{-2}}$ s$^{{-1}}$)")

# difference colorbar
levels_diff = np.linspace(difflim[species][0], difflim[species][1])
m_diff = plt.cm.ScalarMappable(cmap=cmap_diff)
m_diff.set_array(levels_diff)
m_diff.set_clim(difflim[species])

color_bar2 = fig.colorbar(m_diff, orientation="vertical", extend="both", ax=ax[1, :], shrink=0.9, pad=0.005)
color_bar2.set_label(
    f"{print_type} difference \n{species_print[species]} flux {date}\n(mol m$^{{-2}}$ s$^{{-1}}$)"
)


# ### 7. Compare region definitions
#
# THIS NEEDS UPDATING TO COMPARE MULTIPLE MODELS - CURRENTLY THESE NEED TO BE INPUT MANUALLY


species = "ch4"  # select the species you want to plot
plot_region = "UK"  # choose the region mask to plot


# Read in emissions output

intem = xr.open_dataset(
    os.path.join(data_dir, f"InTEM_NAME_EUROPE_{intem_species[species]}_{period[m][species]}.nc")
)
# rhime = xr.open_dataset(os.path.join(data_dir,f'RHIME_NAME_EUROPE_{species}_{period}.nc'))
# empa = xr.open_dataset(os.path.join(data_dir,f'{empa_name}_FLEXPART_EUROPE_{species}_{period}.nc'))

try:
    intem_r0 = np.where(intem["countrynames"].values == plot_region)[0][0]
except:
    print(f"No region in InTEM called {plot_region}")


ax_limits = [-12, 25, 40, 65]  # min_lon, max_lon, min_lat, max_lat

fig, ax = plt.subplots(
    2, 3, figsize=(12, 6), constrained_layout=True, subplot_kw={"projection": cartopy.crs.PlateCarree()}
)

for i in range(2):
    for j in range(3):
        ax[i, j].add_feature(cartopy.feature.BORDERS, linestyle=":", edgecolor="black", linewidth=1.0)
        ax[i, j].coastlines(resolution="50m", color="black", linewidth=1.0)
        ax[i, j].set_extent(ax_limits)

# plot intem
ax[0, 0].pcolormesh(
    intem.lon.values,
    intem.lat.values,
    intem["region_definitions"][:, :, intem_r0],
    cmap="Blues",
    vmin=0,
    vmax=2,
)
ax[0, 0].set_title(f"InTEM {plot_region} mask")

# plot difference
# ax[1,0].pcolormesh(intem.lon.values,intem.lat.values,
#                   intem['region_definitions'][:,:,intem_r0]-rhime['region_definitions'][:,:,rhime_r0],
#                   cmap='bwr')
ax[1, 0].set_title(f"InTEM {plot_region} - RHIME {plot_region}")
