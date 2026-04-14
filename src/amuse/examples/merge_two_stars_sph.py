"""
Module to merge two stars with the use of an SPH code.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from amuse.datamodel import Particle
from amuse.units import nbody_system, units, constants
from amuse.ext.star_to_sph import convert_stellar_model_to_sph
from amuse.community.gadget2 import Gadget2
from amuse.community.evtwin import Evtwin


# #BOOKLISTSTART1# #
def return_evolved_star_hydro(mass, time, number_of_sph_particles):
    star = Particle(mass=mass)
    stellar = Evtwin()
    star = stellar.particles.add_particle(star)
    stellar.evolve_model(time)
    number_of_sph_particles = number_of_sph_particles * int(mass.value_in(units.MSun))
    star_in_sph = convert_stellar_model_to_sph(
        star, number_of_sph_particles
    ).gas_particles
    stellar.stop()
    return star_in_sph


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
def merge_two_stars_sph(
    mass_primary, mass_secondary, time_collision, number_of_sph_particles
):
    primary_in_sph = return_evolved_star_hydro(
        mass_primary, time_collision, number_of_sph_particles
    )
    primary_in_sph = relax_sph_realization(primary_in_sph)
    secondary_in_sph = return_evolved_star_hydro(
        mass_secondary, time_collision, number_of_sph_particles
    )
    secondary_in_sph = relax_sph_realization(secondary_in_sph)
    radius = primary_in_sph.x.max() + secondary_in_sph.x.max()
    total_mass = primary_in_sph.mass.sum() + secondary_in_sph.mass.sum()
    secondary_in_sph.x += 0.8 * radius
    secondary_in_sph.y += 0.6 * radius
    secondary_in_sph.vx -= (constants.G * total_mass / radius).sqrt()

    converter = nbody_system.nbody_to_si(mass_primary, 1.0 | units.au)
    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(primary_in_sph)
    hydro.gas_particles.add_particles(secondary_in_sph)
    hydro.evolve_model(2.0 | units.hour)
    hydro.gas_particles.new_channel_to(primary_in_sph).copy()
    hydro.gas_particles.new_channel_to(secondary_in_sph).copy()
    hydro.stop()
    return primary_in_sph, secondary_in_sph


# #BOOKLISTSTOP2# #


# #BOOKLISTSTART3# #
def relax_sph_realization(sph_star):
    dynamical_timescale = sph_star.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1 | units.RSun)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.gas_particles.add_particles(sph_star)

    to_framework = hydro.gas_particles.new_channel_to(sph_star)

    ts_factor = 2.5
    time_end = ts_factor * sph_star.dynamical_timescale(mass_fraction=0.9)
    n_steps = ts_factor * 100
    velocity_damp_factor = 1.0 - (ts_factor * 2 * np.pi) / n_steps
    time_step = time_end / float(n_steps)
    time = 0 | units.day
    while time < time_end:
        time += time_step
        hydro.evolve_model(time)
        hydro.gas_particles.velocity = (
            velocity_damp_factor * hydro.gas_particles.velocity
        )
    to_framework.copy()
    hydro.stop()
    return sph_star


# #BOOKLISTSTOP3# #


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-M",
        type=units.MSun,
        dest="mass_primary",
        default=10 | units.MSun,
        help="Mass of the primary star",
    )
    result.add_argument(
        "-m",
        type=units.MSun,
        dest="mass_secondary",
        default=1 | units.MSun,
        help="Mass of the secondary star",
    )
    result.add_argument(
        "-N",
        dest="number_of_sph_particles",
        type=int,
        default=100,
        help="Number of sph particles per MSun",
    )
    result.add_argument(
        "-t",
        type=units.Myr,
        dest="time_collision",
        default=0.01 | units.Myr,
        help="end time of the simulation",
    )
    return result


def main(args=sys.argv[1:]):
    arguments = new_argument_parser().parse_args(args)
    p, s = merge_two_stars_sph(**arguments.__dict__)
    plt.scatter(p.x.value_in(units.RSun), p.y.value_in(units.RSun), c="b")
    plt.scatter(s.x.value_in(units.RSun), s.y.value_in(units.RSun), c="r")
    plt.show()


if __name__ == "__main__":
    main()
