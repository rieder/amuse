"""
Example for animating a gravity simulation.
"""

import argparse
import numpy as np

# #BOOKLISTSTART# #
import matplotlib.pyplot as plt
from matplotlib import animation
from amuse.io import read_set_from_file
from amuse.units import units


def animate(x, y):

    def update(i):
        data = np.stack([x[i], y[i]]).T
        scat.set_offsets(data)

        return (scat,)

    number_of_snaps = len(x)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim((-1.2, 1.2))
    ax.set_ylim((-1.2, 1.2))
    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")

    colormap = ["#FFFF00", "wheat", "deepskyblue"]
    size = [40, 20, 20]
    edgecolor = ["orange", "wheat", "deepskyblue"]

    ax.set_facecolor("black")
    scat = ax.scatter(x[0], y[0], c=colormap, s=size, edgecolor=edgecolor)
    anim = animation.FuncAnimation(fig, update, frames=number_of_snaps, interval=30)
    anim.save("anim_gravity.mp4", dpi=150, writer=animation.FFMpegWriter(fps=25))


# #BOOKLISTSTOP# #


def anim_gravity(filename="sun_venus_earth.amuse", **kwargs):
    particles = read_set_from_file(filename)

    x = []
    y = []
    for si in particles.history:
        x.append(si.x.value_in(units.au))
        y.append(si.y.value_in(units.au))

    animate(x, y)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--filename",
        type=str,
        default="sun_venus_earth.amuse",
        help="input filename",
    )
    return parser


def main(**kwargs):
    anim_gravity(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
