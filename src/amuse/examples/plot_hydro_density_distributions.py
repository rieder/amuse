"""
Plots different initial particle distributions for SPH codes
"""

import argparse
from matplotlib import pyplot as plt
import numpy as np

from amuse.units import units, nbody_system
from amuse.datamodel import Particle
from amuse.ic.molecular_cloud import molecular_cloud
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.star_to_sph import convert_stellar_model_to_sph
from amuse.plot import sph_particles_plot
from amuse.community.evtwin import Evtwin
from amuse.community.fi import Fi
from amuse.ext.halogen_model import new_halogen_model


# #BOOKLISTSTART1# #
def stellar_model(
    number_of_particles=10000, mass_star=2.0 | units.MSun, time=0.0 | units.Myr
):
    star = Particle(mass=mass_star)
    stellar_evolution = Evtwin()
    se_star = stellar_evolution.particles.add_particle(star)
    print(f"Evolving {star.mass} to t= {time.in_(units.Myr)}")
    stellar_evolution.evolve_model(time)
    print(f"Stellar type: {se_star.stellar_type}")
    print("Creating SPH particles from the (1D) stellar evolution model")
    sph_particles = convert_stellar_model_to_sph(
        se_star, number_of_particles
    ).gas_particles
    stellar_evolution.stop()
    return sph_particles


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
def plot_zams_stellar_model(number_of_particles=10000, mass_star=2.0 | units.MSun):
    length_unit = units.RSun
    sph_particles = stellar_model(number_of_particles, mass_star)
    figure = plt.figure(figsize=(6, 6))
    ax = figure.add_subplot(1, 1, 1)
    sph_particles_plot(
        sph_particles,
        min_size=500,
        max_size=500,
        alpha=0.01,
        view=(-2, 2, -2, 2) | length_unit,
    )
    ax.set_facecolor("white")
    ax.set_xlabel(f"x [{length_unit}]")
    ax.set_ylabel(f"y [{length_unit}]")

    save_file = "stellar_2MSunZAMS_projected.pdf"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}\n")


# #BOOKLISTSTOP2# #


# #BOOKLISTSTART3# #
def gmc_model(
    number_of_particles=10000,
    mass_cloud=10000.0 | units.MSun,
    radius_cloud=10.0 | units.parsec,
):
    converter = nbody_system.nbody_to_si(mass_cloud, radius_cloud)
    sph_particles = molecular_cloud(
        target_number_of_particles=number_of_particles,
        convert_nbody=converter,
    ).result
    sph = Fi(converter)
    sph.gas_particles.add_particles(sph_particles)
    sph.evolve_model(1 | units.day)
    channel = sph.gas_particles.new_channel_to(sph_particles)
    channel.copy()
    sph.stop()
    return sph_particles


# #BOOKLISTSTOP3# #
def plot_circumstellar_disk_model(
    number_of_particles,
    mass_disk=1.0 | units.MSun,
    radius_disk=100.0 | units.au,
    length_unit=units.au,
):
    converter = nbody_system.nbody_to_si(mass_disk, radius_disk)

    sph_particles = ProtoPlanetaryDisk(
        number_of_particles,
        convert_nbody=converter,
        densitypower=1.5,
        Rmin=0.1,
        Rmax=1,
        q_out=1.0,
        discfraction=0.1,
    ).result
    sph_particles.h_smooth = 0.1 | units.au
    sph_particles.rotate(0.0, np.pi / 6, np.pi / 3)

    figure = plt.figure(figsize=(6, 6))
    ax = figure.add_subplot(111)
    sph_particles_plot(
        sph_particles,
        min_size=500,
        max_size=500,
        alpha=0.01,
        view=(-100, 100, -100, 100) | length_unit,
    )
    ax.set_facecolor("white")
    ax.set_xlabel(f"x [{length_unit}]")
    ax.set_ylabel(f"y [{length_unit}]")

    save_file = "fig_circum_stellar_disk_model.pdf"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}\n")


def plot_gmc_model(
    number_of_particles=10000,
    mass_cloud=10000.0 | units.MSun,
    radius_cloud=10.0 | units.parsec,
    length_unit=units.pc,
):
    sph_particles = gmc_model(number_of_particles, mass_cloud, radius_cloud)
    figure = plt.figure(figsize=(6, 6))
    ax = figure.add_subplot(111)
    sph_particles_plot(
        sph_particles,
        min_size=500,
        max_size=500,
        alpha=0.01,
        view=(-20, 20, -20, 20) | length_unit,
    )
    ax.set_facecolor("white")
    ax.set_xlabel(f"x [{length_unit}]")
    ax.set_ylabel(f"y [{length_unit}]")

    file = "molecular_cloud_projected.pdf"
    plt.savefig(file)
    print("Saved figure in file", file)


def sph_galaxy_plot(
    particles,
    u_range=None,
    min_size=100,
    max_size=10000,
    alpha=0.1,
    gd_particles=None,
    width=None,
    view=None,
    length_unit=units.kpc,
    speed_unit=units.kms,
):
    x = particles.x.value_in(length_unit)
    y = particles.y.value_in(length_unit)
    c = (particles.velocity.lengths().value_in(speed_unit)) ** 2
    plt.scatter(x, y, s=500, c=c, edgecolors="none", alpha=alpha, cmap="jet_r")
    plt.xlabel(f"x [{length_unit}]")
    plt.ylabel(f"y [{length_unit}]")


def plot_disk_galaxy_model(
    number_of_particles=10000,
    mass_galaxy=1e10 | units.MSun,
    radius_galaxy=30 | units.kpc,
    x_label="x [length]",
    y_label="y [length]",
):
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_aspect("equal")
    plt.xlim(-100, 100)
    plt.ylim(-100, 100)

    converter = nbody_system.nbody_to_si(mass_galaxy, radius_galaxy)
    sph_particles = new_halogen_model(
        number_of_particles,
        convert_nbody=converter,
        alpha=1.0,
        beta=5.0,
        gamma=0.5,
    )
    # sph_particles.u = sph_particles.velocity.lengths()**2
    sph_particles.u = 1 | units.kms**2
    sph_particles[:10].u = sph_particles[:10].velocity.lengths() ** 2
    sph_particles.h_smooth = 0.1 | units.kpc
    sph_galaxy_plot(
        sph_particles,
        min_size=500,
        max_size=500,
        alpha=0.01,
        view=(-100, 100, -100, 100) | units.kpc,
    )

    save_file = "disk_galaxy_model.pdf"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}")


def plot_hydro_density_distributions(
    number_of_particles=10000,
    mass_star=2.0 | units.MSun,
    radius_cloud=10.0 | units.parsec,
    mass_cloud=10000.0 | units.MSun,
    radius_disk=100.0 | units.au,
    mass_disk=1.0 | units.MSun,
    radius_galaxy=30 | units.kpc,
    mass_galaxy=1e10 | units.MSun,
):
    print("plotting")
    plot_zams_stellar_model(number_of_particles, mass_star=mass_star)
    plot_gmc_model(
        number_of_particles, mass_cloud=mass_cloud, radius_cloud=radius_cloud
    )
    plot_circumstellar_disk_model(
        number_of_particles, mass_disk=mass_disk, radius_disk=radius_disk
    )
    plot_disk_galaxy_model(
        number_of_particles, mass_galaxy=mass_galaxy, radius_galaxy=radius_galaxy
    )
    print("")


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-n",
        "--number_of_particles",
        type=int,
        default=10000,
        help="Number of particles",
    )
    parser.add_argument(
        "-ms",
        "--mass_star",
        type=units.MSun,
        default=10000.0 | units.MSun,
        help="Mass of the star",
    )
    parser.add_argument(
        "-rc",
        "--radius_cloud",
        type=units.parsec,
        default=10.0 | units.parsec,
        help="Radius of the cloud",
    )
    parser.add_argument(
        "-mc",
        "--mass_cloud",
        type=units.MSun,
        default=10000.0 | units.MSun,
        help="Mass of the cloud",
    )
    parser.add_argument(
        "-rd",
        "--radius_disk",
        type=units.au,
        default=100.0 | units.au,
        help="Radius of the disk",
    )
    parser.add_argument(
        "-md",
        "--mass_disk",
        type=units.MSun,
        default=1.0 | units.MSun,
        help="Mass of the disk",
    )
    parser.add_argument(
        "-rg",
        "--radius_galaxy",
        type=units.kpc,
        default=30.0 | units.kpc,
        help="Radius of the galaxy",
    )
    parser.add_argument(
        "-mg",
        "--mass_galaxy",
        type=units.MSun,
        default=1e10 | units.MSun,
        help="Mass of the galaxy",
    )
    return parser


def main(**kwargs):
    plot_hydro_density_distributions(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
