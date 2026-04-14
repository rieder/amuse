"""
Salpeter mass function for 1000 stars between 1 MSun and 100 MSun, overplotted
with an analytic power law with a slope of -2.35.
"""

import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt

from amuse.ic.salpeter import new_powerlaw_mass_distribution
from amuse.ic.brokenimf import new_kroupa_mass_distribution
from amuse.ic.brokenimf import new_miller_scalo_mass_distribution
from amuse.units import units


# #BOOKLISTSTART# #
def plot_mass_function(masses, ximf, color, label, n_bins=51):
    "Plots a histogram of the masses provided, with a power law added if ximf > 0"
    mass_min = masses.min()
    mass_max = masses.max()
    log_min = np.log10(mass_min.value_in(units.MSun))
    log_max = np.log10(mass_max.value_in(units.MSun))
    bins = np.logspace(log_min, log_max, n_bins)
    print(bins)
    bin_number, bin_edges = np.histogram(masses.value_in(units.MSun), bins=bins)
    y = bin_number / (bin_edges[1:] - bin_edges[:-1])
    x = (bin_edges[1:] + bin_edges[:-1]) / 2.0
    for yi in y:
        yi = max(yi, 1.0e-10)

    plt.scatter(x, y, s=50, c=color, lw=0, label=label)

    if ximf < 0:
        c = (
            (mass_max.value_in(units.MSun) ** (ximf + 1))
            - (mass_min.value_in(units.MSun) ** (ximf + 1))
        ) / (ximf + 1)
        plt.plot(x, len(masses) / c * (x**ximf), c=color)


# #BOOKLISTSTOP# #


def plot_power_law_mass_functions(
    number_of_stars=10000,
    mass_min=0.1 | units.MSun,
    mass_max=100 | units.MSun,
    ximf=2.35,
    seed=31415,
):
    """
    Plots Miller-Scalo, Kroupa and single-power-law (Salpeter) mass functions.
    """
    np.random.seed(seed)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"M [M$_\odot$]")
    ax.set_ylabel("N")
    ax.set_xlim([1e-1, 1e2])
    ax.set_ylim([1e-2, 1e5])
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    masses = new_miller_scalo_mass_distribution(number_of_stars, mass_min, mass_max)
    plot_mass_function(masses, +1, color[2], "Miller-Scalo")

    masses = new_kroupa_mass_distribution(number_of_stars, mass_min, mass_max)
    plot_mass_function(masses, +1, color[1], "Kroupa")

    masses = new_powerlaw_mass_distribution(number_of_stars, mass_min, mass_max, ximf)
    plot_mass_function(masses, ximf, color[0], "Salpeter")

    plt.legend()
    save_file = "salpeter.pdf"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}\n")


def new_argument_parser(args):
    """
    Parse provided arguments
    """
    result = argparse.ArgumentParser(
        args,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-N", "--number_of_stars", type=int, default=10000, help="number of stars"
    )
    result.add_argument(
        "-m",
        "--mass_min",
        type=units.MSun,
        default=0.1 | units.MSun,
        help="minimum mass of the mass function",
    )
    result.add_argument(
        "-M",
        "--mass_max",
        type=units.MSun,
        default=100 | units.MSun,
        help="maximum mass of the mass function",
    )
    result.add_argument(
        "-x", "--ximf", type=float, default=-2.35, help="mass function slope"
    )
    result.add_argument("-s", "--seed", type=int, default=31415, help="random seed")
    return result


def salpeter(**kwargs):
    plot_power_law_mass_functions(**kwargs)


def main(**kwargs):
    salpeter(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser(sys.argv[1:]).parse_args()
    main(**arguments.__dict__)
