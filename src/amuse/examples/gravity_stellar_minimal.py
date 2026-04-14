"""
AMUSE example combining gravity and stellar evolution.
"""

import sys
import argparse
from amuse.examples import shell_is_interactive
from amuse.units import units, nbody_system
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.kingmodel import new_king_model
from amuse.io import write_set_to_file
from amuse.community.sse import Sse
from amuse.community.ph4 import Ph4


def gravity_stellar_minimal(
    number_of_stars,
    time_end,
    W0,
    radius_virial,
    mass_min,
    mass_max,
    metallicity,
    filename="gravity_stellar_minimal.amuse",
):
    """
    Sets up a simulation that combines gravity and stellar evolution, and runs
    this until 'time_end', writing the result to a file.
    """
    # #BOOKLISTSTART1# #
    masses = new_salpeter_mass_distribution(number_of_stars, mass_min, mass_max)
    mass_total_init = masses.sum()
    converter = nbody_system.nbody_to_si(mass_total_init, radius_virial)
    bodies = new_king_model(number_of_stars, W0, convert_nbody=converter)
    bodies.mass = masses
    bodies.scale_to_standard(convert_nbody=converter)

    stellar = Sse()
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particle(bodies)
    bodies.radius = stellar.particles.radius

    gravity = Ph4(converter)
    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)

    channel_from_stellar_to_gravity = stellar.particles.new_channel_to(
        gravity.particles
    )
    channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)
    # #BOOKLISTSTOP1# #

    # #BOOKLISTSTART2# #
    energy_total_init = gravity.kinetic_energy + gravity.potential_energy
    energy_difference_gr = 0 | energy_total_init.unit

    time = 0.0 | time_end.unit
    while time < time_end:
        time_step = min(stellar.particles.time_step.amin(), time_end - time)

        stellar.evolve_model(time + time_step / 2)
        channel_from_stellar_to_gravity.copy()

        energy_total_gr = gravity.kinetic_energy + gravity.potential_energy
        gravity.evolve_model(time + time_step)
        energy_difference_gr += (
            gravity.kinetic_energy + gravity.potential_energy - energy_total_gr
        )

        stellar.evolve_model(time + time_step)
        channel_from_stellar_to_gravity.copy()
        channel_from_gravity_to_framework.copy()

        time += time_step
        write_set_to_file(bodies, filename, append_to_file=True)
    # #BOOKLISTSTOP2# #

    energy_kinetic = gravity.kinetic_energy
    energy_potential = gravity.potential_energy
    energy_total = energy_kinetic + energy_potential
    energy_difference = energy_total_init - energy_total
    mass_total = bodies.mass.sum()
    print(
        f"T= {time.in_(units.Myr)} "
        f"M= {mass_total.in_(units.MSun)} (dM[SE]={mass_total / mass_total_init}) "
        f"E= {energy_total.in_(units.erg)} Q= {energy_kinetic / energy_potential} "
        f"dE/E= {(energy_total_init - energy_total) / energy_total} "
        f"(dE[gr]/E= {energy_difference_gr / energy_total}, "
        f"dE[se]/E= "
        f"{(energy_total_init - energy_total - energy_difference_gr) / energy_total})"
    )
    energy_total_init -= energy_difference

    gravity.stop()
    stellar.stop()


def gravity_stellar_minimal_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-N",
        "--number_of_stars",
        dest="number_of_stars",
        type=int,
        default=100,
        help="number of stars",
    )
    result.add_argument(
        "-M",
        "--mass_max",
        type=units.MSun,
        dest="mass_max",
        default=100 | units.MSun,
        help="maximal stellar mass",
    )
    result.add_argument(
        "-m",
        "--mass_min",
        type=units.MSun,
        dest="mass_min",
        default=0.1 | units.MSun,
        help="minimal stellar mass",
    )
    result.add_argument(
        "-R",
        "--radius_virial",
        type=units.parsec,
        dest="radius_virial",
        default=1.0 | units.parsec,
        help="cluster virial radius",
    )
    result.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        dest="time_end",
        default=1.0 | units.Myr,
        help="end time of the simulation",
    )
    result.add_argument(
        "-W",
        "--king_w0",
        dest="W0",
        type=float,
        default=7.0,
        help="structure parameter for King model",
    )
    result.add_argument(
        "-z",
        "--metallicity",
        dest="metallicity",
        type=float,
        default=0.02,
        help="metallicity",
    )
    result.add_argument(
        "-f",
        "--filename",
        dest="filename",
        type=str,
        default="gravity_stellar_minimal.amuse",
        help="name of the file containing the results",
    )
    return result


def main():
    if shell_is_interactive():
        gravity_stellar_minimal({})
    else:
        arguments = gravity_stellar_minimal_argument_parser().parse_args(sys.argv[1:])
        gravity_stellar_minimal(**arguments.__dict__)


if __name__ == "__main__":
    main()
