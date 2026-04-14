"""
Example AMUSE script to generate a Plummer sphere and plot the results.
"""

import argparse

# #BOOKLISTSTART# #
import matplotlib.pyplot as plt
from amuse.units.nbody_system import length
from amuse.ic.plummer import new_plummer_model


def plot_plummer(number_of_particles=1000):
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(1, 1, 1)
    bodies = new_plummer_model(number_of_particles)
    ax.scatter(bodies.x.value_in(length), bodies.y.value_in(length))
    ax.set_xlim((-1, 1))
    ax.set_ylim((-1, 1))
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    plt.savefig("plummer.pdf")
    print("Saved figure in file plummer.pdf")


# #BOOKLISTSTOP# #


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-N",
        "--number_of_particles",
        type=int,
        default=1000,
        help="number of stars",
    )
    return result


def main(**kwargs):
    plot_plummer(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
