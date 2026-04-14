"""
   Minimal routine for running a hydrodynamics code
"""

import sys
import argparse
from amuse.units import nbody_system, units
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.gadget2 import Gadget2
from amuse.io import write_set_to_file


# #BOOKLISTSTART1# #
def hydro_minimal(
    number_of_particles=100,
    total_mass=1 | units.MSun,
    virial_radius=1 | units.RSun,
    time_end=6 | units.hour,
):
    converter = nbody_system.nbody_to_si(total_mass, virial_radius)
    gas = new_plummer_gas_model(number_of_particles, convert_nbody=converter)

    hydro = Gadget2(converter)
    hydro.gas_particles.add_particles(gas)
    energy_total_init = (
        hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    )
    hydro.evolve_model(time_end)
    write_set_to_file(hydro.particles, "hydro.amuse")

    energy_kinetic = hydro.kinetic_energy
    energy_potential = hydro.potential_energy
    energy_thermal = hydro.thermal_energy
    energy_total = energy_kinetic + energy_potential + energy_thermal
    virial_ratio = (energy_kinetic + energy_thermal) / energy_potential
    energy_difference = (energy_total_init - energy_total) / energy_total
    com = hydro.gas_particles.center_of_mass()
    print(
        f"T= {hydro.get_time()} M= {hydro.gas_particles.mass.sum()} "
        f"E= {energy_total} Q= {virial_ratio} "
        f"dE= {energy_difference} CoM= {com.in_(units.RSun)}"
    )

    hydro.stop()
# #BOOKLISTSTOP1# #


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-N",
        "--number_of_particles",
        type=int,
        default=100,
        help="number of gas particles",
    )
    result.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=6 | units.hour,
        help="end time of the simulation",
    )
    result.add_argument(
        "-M",
        "--total_mass",
        type=units.MSun,
        default=1 | units.MSun,
        help="Mass of the cloud",
    )
    result.add_argument(
        "-R",
        "--virial_radius",
        type=units.RSun,
        default=1 | units.RSun,
        help="Radius of the cloud",
    )
    return result


def main(**kwargs):
    hydro_minimal(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
