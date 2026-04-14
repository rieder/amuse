"""
Minimal routine to plot N-body simulation data.
"""

import sys
import argparse
from amuse.units import units

# #BOOKLISTSTART# #
import matplotlib.pyplot as plt
from amuse.io import read_set_from_file


def plot_gravity(
    filename="sun_venus_earth.amuse",
    length_unit=units.au,
):
    particles = read_set_from_file(filename)

    figure = plt.figure(figsize=(8, 8))
    ax = figure.add_subplot(1, 1, 1)

    colormap = ["yellow", "green", "blue"]  # specific to a 3-body plot
    size = [40, 20, 20]
    edgecolor = ["orange", "green", "blue"]

    for i, snap in enumerate(particles.history):
        ax.scatter(
            snap.x.value_in(length_unit),
            snap.y.value_in(length_unit),
            c=colormap,
            s=size,
            edgecolor=edgecolor,
        )
    ax.set_xlabel(f"x [{length_unit}]")
    ax.set_ylabel(f"y [{length_unit}]")

    save_file = "plot_gravity.png"
    plt.savefig(save_file)
    print(f"\nSaved figure in file {save_file}\n")


# #BOOKLISTSTOP# #


def new_argument_parser():
    result = argparse.ArgumentParser()
    result.add_argument(
        "-i", "--filename", type=str, default="sun_venus_earth.amuse", help="input file"
    )
    return result


def main(**kwargs):
    plot_gravity(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
