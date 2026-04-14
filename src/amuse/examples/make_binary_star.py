"""
Example for creating a binary star particle
"""

import sys
import argparse
import numpy as np
from amuse.datamodel import Particles, Particle
from amuse.units import constants, units, nbody_system
from amuse.units.quantities import as_vector_quantity
from amuse.community.kepler import Kepler


def semi_major_axis_to_orbital_period(a, mass_total):
    """
    Takes a semi-major axis and total mass, returns the orbital period
    """
    return 2 * np.pi * (a**3 / (constants.G * mass_total)).sqrt()


# #BOOKLISTSTART2# #
def new_kepler(converter):
    """Starts and returns a new Kepler integrator with given converter"""
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.set_longitudinal_unit_vector(1.0, 0.0, 0.0)
    kepler.set_transverse_unit_vector(0.0, 1.0, 0)
    return kepler


# #BOOKLISTSTOP2# #


# #BOOKLISTSTART1# #
def new_binary_orbit(stars, kepler):
    """
    Takes a Particles object of length 2 and returns a binary particle from the
    stars, using the values from `kepler`
    """
    if len(stars) != 2:
        raise ValueError("'stars' must consist of exactly two particles")
    rel_position = as_vector_quantity(kepler.get_separation_vector())
    rel_velocity = as_vector_quantity(kepler.get_velocity_vector())

    mu = stars[0].mass / stars.mass.sum()
    binary = Particle()
    binary.child1 = stars[0]
    binary.child2 = stars[1]
    binary.child1.position = mu * rel_position
    binary.child2.position = -(1 - mu) * rel_position
    binary.child1.velocity = -(1 - mu) * rel_velocity
    binary.child2.velocity = mu * rel_velocity
    print(binary.child1.position.in_(units.au))
    print(binary.child2.position.in_(units.au))
    return binary


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART3# #
def calculate_orbital_elements(bi, kepler):
    """
    Takes a binary particle and a kepler code, returns the semi-major axis and
    eccentricity
    """
    comp1 = bi.child1
    comp2 = bi.child2
    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kepler.initialize_from_dyn(
        mass,
        pos[0],
        pos[1],
        pos[2],
        vel[0],
        vel[1],
        vel[2],
    )
    a, e = kepler.get_elements()
    return a, e


# #BOOKLISTSTOP3# #


def orbital_elements_from_amuse(binary):
    from amuse.ext.orbital_elements import orbital_elements

    stars = Particles()
    stars.add_particle(binary.child1)
    stars.add_particle(binary.child2)
    a, e = orbital_elements(stars)[2:4]
    return a, e


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-M",
        "--mass_1",
        type=units.MSun,
        default=2.0 | units.MSun,
        help="mass of the primary star",
    )
    result.add_argument(
        "-m",
        "--mass_2",
        type=units.MSun,
        default=1.0 | units.MSun,
        help="mass of the secondary star",
    )
    result.add_argument(
        "-a",
        "--semi_major_axis",
        type=units.au,
        default=5.2 | units.au,
        help="separation",
    )
    result.add_argument(
        "-e", "--eccentricity", type=float, default=0.6, help="binary eccentricity"
    )
    return result


def make_binary_star(
    mass_1=2.0 | units.MSun,
    mass_2=1.0 | units.MSun,
    semi_major_axis=5.2 | units.au,
    eccentricity=0.6,
    kepler=None,
):
    stars = Particles(2)
    stars[0].mass = mass_1
    stars[1].mass = mass_2
    converter = nbody_system.nbody_to_si(
        stars.mass.sum(),
        semi_major_axis,
    )
    kepler = new_kepler(converter)
    kepler.initialize_from_elements(
        stars.total_mass(),
        semi_major_axis,
        eccentricity,
    )
    binary = new_binary_orbit(stars, kepler)
    a, e = orbital_elements_from_amuse(binary)
    print(f"Orbital parameters: {a.in_(units.au)} {e:4.2f}")
    kepler.stop()


def main(**kwargs):
    make_binary_star(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args(sys.argv[1:])
    main(**arguments.__dict__)
