"""
Plot the core temperature and density of a star
"""

import sys
import argparse
import pickle
import matplotlib.pyplot as plt

from amuse.units import units

from .plot_stellar_evolution_track import get_color_based_on_stellar_type

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
HeWhiteDwarf = 10 | units.stellar_type


def plot_as_line_per_stellar_type(rhoc, Tc, stp):
    # palette = sns.color_palette(n_colors=25)
    # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']

    ilast = 0
    for inext in range(len(rhoc)):
        if stp[ilast] != stp[inext]:
            # color = palette[stp[ilast].value_in(units.stellar_type)]
            color = get_color_based_on_stellar_type(stp[ilast])
            plt.plot(
                rhoc[ilast:inext].value_in(units.g / units.cm**3),
                Tc[ilast:inext].value_in(units.K),
                c=color,
                lw=4,
            )
            ilast = inext
    inext = -1
    # color = palette[stp[ilast].value_in(units.stellar_type)]
    color = get_color_based_on_stellar_type(stp[ilast])
    plt.plot(
        rhoc[ilast:inext].value_in(units.g / units.cm**3),
        Tc[ilast:inext].value_in(units.K),
        c=color,
        lw=4,
    )


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-f", dest="filenames", nargs="+", default="", help="Stellar data files"
    )
    return result


def plot_core_temperature_density(filenames=""):
    # mass_list = [1, 10, 100] | units.MSun
    if filenames == "":
        filenames = [
            "stellar_Mesa_r2208_M1.0MSun_core_temperature_and_density.pkl",
            "stellar_Mesa_r2208_M10.0MSun_core_temperature_and_density.pkl",
            "stellar_Mesa_r2208_M100.0MSun_core_temperature_and_density.pkl",
        ]

    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_xscale("log")
    ax.set_yscale("log")

    stps = []
    for fi in filenames:
        time, rhoc, Tc, stp = pickle.load(open(fi, "rb"))
        for stpi in stp:
            if stpi not in stps:
                stps.append(stpi)
        plot_as_line_per_stellar_type(rhoc, Tc, stp)

    # palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
    # palette = sns.color_palette(n_colors=25)
    for stpi in stps:
        color = get_color_based_on_stellar_type(stpi)
        # color = palette[stpi.value_in(units.stellar_type)]
        ax.scatter([1], [1], c=color, label=stpi)

    fontsize = 12
    ax.legend(loc="upper left", fontsize=8)
    ax.text(20.0, 1.5e7, r"$1\,M_\odot$", fontsize=fontsize)
    ax.text(1.0e2, 7.0e7, r"$10\,M_\odot$", fontsize=fontsize)
    ax.text(4.0e4, 1.6e9, r"$100\,M_\odot$", fontsize=fontsize)

    ax.set_xlabel("core density [g/cm$^3$]")
    ax.set_ylabel("core temperature [K]")
    ax.set_xlim((1.0, 1.0e8))
    ax.set_ylim((3.0e6, 1.0e10))

    save_file = "plot_core_temperature_density.pdf"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}")


def main(**kwargs):
    plot_core_temperature_density(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    plot_core_temperature_density(**arguments.__dict__)
