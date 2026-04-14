"""
Initialize two stars to a certain age and merge them using MMAMS
"""

import sys
import argparse
import matplotlib.pyplot as plt
from amuse.datamodel import Particle
from amuse.units import units
from amuse.community.mesa import Mesa
from amuse.community.mmams import Mmams


# #BOOKLISTSTART# #
def merge_two_stars(
    mass_primary=5 | units.MSun,
    mass_secondary=3 | units.MSun,
    time_collision=1.0 | units.Myr,
):
    """
    Merges two stars using the Make Me a Massive Star code and returns the
    resulting density profile.
    """
    primary = Particle(mass=mass_primary)
    secondary = Particle(mass=mass_secondary)

    stellar = Mesa()
    primary = stellar.particles.add_particle(primary)
    secondary = stellar.particles.add_particle(secondary)

    stellar.evolve_model(time_collision)

    stellar.merge_colliding(
        primary.copy(), secondary.copy(), Mmams, return_merge_products=["se"]
    )
    radius = stellar.particles[0].get_radius_profile()
    rho = stellar.particles[0].get_density_profile()
    stellar.stop()
    plot_density_profile(radius, rho)


# #BOOKLISTSTOP# #


def plot_density_profile(radius, rho):
    """
    Plots density against radius.
    """
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.plot(radius.value_in(units.RSun), rho.value_in(units.g / units.cm**3))
    ax.set_xlabel(r"$R$ [$R_\odot$]")
    ax.set_ylabel("density [$g/cm^3$]")
    ax.set_yscale("log", nonpositive="clip")
    plt.savefig("merge_two_stars.pdf")
    print("Saved figure in file merge_two_stars.pdf")


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-M",
        "--mass_primary",
        type=units.MSun,
        default=5 | units.MSun,
        help="Mass of the primary star",
    )
    result.add_argument(
        "-m",
        "--mass_secondary",
        type=units.MSun,
        default=3 | units.MSun,
        help="Mass of the secondary star",
    )
    result.add_argument(
        "-t",
        "--time_collision",
        type=units.Myr,
        default=1.0 | units.Myr,
        help="end time of the simulation",
    )
    return result


def main(**kwargs):
    merge_two_stars(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
