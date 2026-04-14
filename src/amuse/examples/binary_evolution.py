"""
Generates a grid of binaries with different, primary mass, mass ratio
and separation and evolves these over time.
"""

import argparse
import matplotlib.pyplot as plt

from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.seba import Seba


# #BOOKLISTSTART1# #
def create_double_star(mass_primary, mass_secondary, semi_major_axis, eccentricity):
    primary_stars = Particles(mass=mass_primary)
    secondary_stars = Particles(mass=mass_secondary)

    stars = Particles()
    primary_stars = stars.add_particles(primary_stars)
    secondary_stars = stars.add_particles(secondary_stars)

    double_star = Particles(semi_major_axis=semi_major_axis, eccentricity=eccentricity)
    double_star.child1 = list(primary_stars)
    double_star.child2 = list(secondary_stars)
    return double_star, stars


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
def evolve_double_star(
    mass_primary,
    mass_secondary,
    semi_major_axis,
    eccentricity,
    time_end,
    number_of_steps,
):
    double_star, stars = create_double_star(
        mass_primary, mass_secondary, semi_major_axis, eccentricity
    )
    time = 0 | units.Myr
    time_step = time_end / number_of_steps

    code = Seba()
    code.particles.add_particles(stars)
    code.binaries.add_particles(double_star)

    channel_from_code_to_model_for_binaries = code.binaries.new_channel_to(double_star)

    t = [] | units.Myr
    a = [] | units.RSun
    e = []
    while time < time_end:
        time += time_step
        code.evolve_model(time)
        channel_from_code_to_model_for_binaries.copy()
        t.append(time)
        a.append(double_star[0].semi_major_axis)
        e.append(double_star[0].eccentricity)
    code.stop()
    return t, a, e


# #BOOKLISTSTOP2# #


def plot_semi_major_axis_eccentricity(
    t,
    a,
    e,
    time_unit=units.Myr,
    length_unit=units.au,
):
    fig = plt.figure(figsize=(8, 6))
    fta = fig.add_subplot(2, 1, 1)
    fta.plot(t.value_in(time_unit), a.value_in(length_unit))
    plt.ylabel(f"semi major axis ({length_unit})")
    fte = fig.add_subplot(2, 1, 2)
    fte.plot(t.value_in(time_unit), e)
    plt.ylabel("eccentricity")
    plt.xlabel(f"time [{time_unit}]")
    plt.savefig("binary_evolution.png")


def new_argument_parser():
    result = argparse.ArgumentParser()
    result.add_argument(
        "-M",
        "--mass_primary",
        type=units.MSun,
        dest="mass_primary",
        default=12 | units.MSun,
        help="mass of the primary star",
    )
    result.add_argument(
        "-m",
        "--mass_secondary",
        type=units.MSun,
        dest="mass_secondary",
        default=10 | units.MSun,
        help="mass of the secondary star",
    )
    result.add_argument(
        "-T",
        "--time_end",
        type=units.Myr,
        dest="time_end",
        default=25.0 | units.Myr,
        help="end time of the simulation",
    )
    result.add_argument(
        "-a",
        "--semi_major_axis",
        type=units.RSun,
        dest="semi_major_axis",
        default=205 | units.RSun,
        help="orbital separation",
    )
    result.add_argument(
        "-e",
        "--eccentricity",
        dest="eccentricity",
        type=float,
        default=0.0,
        help="orbital eccentricity",
    )
    result.add_argument(
        "-n",
        "--number_of_steps",
        dest="number_of_steps",
        type=int,
        default=100,
        help="number of output steps",
    )
    return result


def binary_evolution(
    mass_primary=12 | units.MSun,
    mass_secondary=10 | units.MSun,
    time_end=25.0 | units.Myr,
    semi_major_axis=205 | units.RSun,
    eccentricity=0.0,
    number_of_steps=100,
):
    t, a, e = evolve_double_star(
        mass_primary=mass_primary,
        mass_secondary=mass_secondary,
        semi_major_axis=semi_major_axis,
        eccentricity=eccentricity,
        time_end=time_end,
        number_of_steps=number_of_steps,
    )
    plot_semi_major_axis_eccentricity(t, a, e)


def main(**kwargs):
    binary_evolution(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
