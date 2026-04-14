"""
Example to plot density distributions for different initial condition generators.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from amuse.units import nbody_system, units
from amuse.datamodel import Particle
from amuse.ic.plummer import new_plummer_model
from amuse.ic.kingmodel import new_king_model
from amuse.ic.fractalcluster import new_fractal_cluster_model
from amuse.ext.halogen_model import new_halogen_model
from amuse.ext.protodisk import ProtoPlanetaryDisk

from amuse_examples.make_oort_cloud import add_comets


def plot_projected_density(
    model,
    figure=None,
    ax=None,
    color_number=0,
    x_axis="x",
    y_axis="z",
    x_unit=nbody_system.length,
    y_unit=nbody_system.length,
    x_lim=(-1, 1),
    y_lim=(-1, 1),
):
    """
    Plots the particles in 'model', returns the figure
    """
    if figure is None:
        figure = plt.figure()
    if ax is None:
        ax = figure.add_subplot(111)
    x = getattr(model, x_axis)
    y = getattr(model, y_axis)
    x_label = f"{x_axis} [{x.unit}]"
    y_label = f"{y_axis} [{y.unit}]"
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_aspect("equal")
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.scatter(
        x.value_in(x_unit),
        y.value_in(y_unit),
        c=color[color_number],
        s=40,
        lw=0,
    )
    return figure


def plot_plummer_model(
    number_of_particles,
    x_axis="x",
    y_axis="y",
):
    """
    Generates a Plummer model and plots its projected density, returns the figure
    """
    model = new_plummer_model(number_of_particles)
    figure = plot_projected_density(
        model,
        color_number=0,
        x_axis=x_axis,
        y_axis=y_axis,
    )
    return figure


def plot_king_model(
    number_of_particles,
    W=9,
    x_axis="x",
    y_axis="y",
):
    """
    Generates a King model and plots its projected density, returns the figure
    """
    model = new_king_model(number_of_particles, W)
    figure = plot_projected_density(
        model,
        color_number=1,
        x_axis=x_axis,
        y_axis=y_axis,
    )
    return figure


def plot_hierarchical_king_model(
    number_of_particles,
    W=9,
    number_of_levels=5,
    x_axis="x",
    y_axis="z",
):
    """
    Generates a hierarchical King model and plots its projected density,
    returns the figure
    """
    model_com = new_king_model(5, W)
    model_com.position *= 0.7
    figure = None
    ax = None
    for i in range(number_of_levels):
        model = new_king_model(int(number_of_particles / number_of_levels), W)
        model.position *= 0.1
        model.position += model_com[i].position
        figure = plot_projected_density(
            model,
            figure=figure,
            ax=ax,
            color_number=2,
            x_axis=x_axis,
            y_axis=y_axis,
        )
        ax = figure.gca()
    return figure


def plot_fractal_model(
    number_of_particles,
    fractal_dimension=1.6,
    x_axis="x",
    y_axis="y",
):
    """
    Generates a fractal model and plots its projected density, returns the figure
    """
    model = new_fractal_cluster_model(
        N=number_of_particles, fractal_dimension=fractal_dimension, random_seed=42
    )
    figure = plot_projected_density(
        model,
        color_number=3,
        x_axis=x_axis,
        y_axis=y_axis,
    )
    return figure


def plot_galaxy_model(
    number_of_particles,
    x_axis="x",
    y_axis="y",
):
    """
    Generates a galaxy model and plots its projected density, returns the figure
    """
    model = new_halogen_model(number_of_particles, alpha=1.0, beta=5.0, gamma=0.5)
    figure = plot_projected_density(model, color_number=4, x_axis=x_axis, y_axis=y_axis)
    return figure


def plot_circumstellar_disk_model(
    number_of_particles,
    x_axis="x",
    y_axis="y",
):
    """
    Generates a circumstellar disk model and plots its projected density,
    returns the figure
    """
    model = ProtoPlanetaryDisk(
        number_of_particles,
        densitypower=1.5,
        Rmin=0.1,
        Rmax=1,
        q_out=1.0,
        discfraction=0.1,
    ).result
    model.rotate(0.0, np.pi / 6, np.pi / 6)

    figure = plot_projected_density(model, color_number=5, x_axis=x_axis, y_axis=y_axis)
    ax = figure.gca()
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.scatter([0], [0], s=200, c=color[8], lw=0)
    return figure


def plot_oort_cloud_around_star_model(number_of_comets, x_axis="x", y_axis="y"):
    """
    Generates an Oort cloud model and plots its projected density, returns the
    figure
    """
    star = Particle(mass=1 | units.MSun)
    star.position = (0, 0, 0) | units.au
    star.velocity = (0, 0, 0) | units.kms
    m_comets = 0.001 * star.mass
    n_comets = number_of_comets
    q_min = 10 | units.au
    a_min = 1000 | units.au
    a_max = 2.0e5 | units.au
    model = add_comets(star, m_comets, n_comets, q_min, a_min, a_max, seed=1)

    figure = plot_projected_density(
        model,
        color_number=5,
        x_axis=x_axis,
        y_axis=y_axis,
        x_lim=(-0.2, 0.2),
        y_lim=(-0.2, 0.2),
        x_unit=units.pc,
        y_unit=units.pc,
    )
    ax = figure.gca()
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax.scatter([0], [0], s=200, c=color[8], lw=0)
    return figure


def main(number_of_particles=10):
    """
    Plots density distributions for a number of initial condition models,
    saves the figures to disk
    """
    np.random.seed(42)
    for model_type in [
        "plummer",
        "king",
        "fractal",
        "galaxy",
        "hierarchical_king",
        # "circumstellar_disk",
        "oort_cloud_around_star",
        # "oligarchic_planetary_system",
    ]:
        plot_function = locals().get(f"plot_{model_type}_model")
        if plot_function is None:
            print(f"Using globals for {model_type}")
            plot_function = globals().get(f"plot_{model_type}_model")

        figure = plot_function(number_of_particles)
        filename = f"{model_type}_model.pdf"
        plt.figure(figure)
        plt.savefig(filename)
        print(f"Saved figure in file {filename}")


def new_argument_parser(args):
    """
    Parses command line arguments
    """
    result = argparse.ArgumentParser(
        args,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-N",
        "--number_of_particles",
        dest="number_of_particles",
        type=int,
        default=1000,
        help="number of particles",
    )
    return result.parse_args()


if __name__ == "__main__":
    arguments = new_argument_parser(sys.argv[1:])
    main(**arguments.__dict__)
