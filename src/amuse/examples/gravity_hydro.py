"""
example code for bridging a gravity solver with a hydrodynamics solver
"""

import argparse
import numpy as np
from amuse.units import nbody_system, units, constants
from amuse.io import write_set_to_file
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.community.fi import Fi
from amuse.community.ph4 import Ph4
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary


# #BOOKLISTSTART1# #
class BaseCode:
    def __init__(self, code, particles, eps=0 | units.RSun):

        self.local_particles = particles
        m = self.local_particles.mass.sum()
        l = self.local_particles.position.length()
        self.converter = nbody_system.nbody_to_si(m, l)
        self.code = code(self.converter)
        self.code.parameters.epsilon_squared = eps**2

    def evolve_model(self, time):
        self.code.evolve_model(time)

    def copy_to_framework(self):
        self.channel_to_framework.copy()

    def get_gravity_at_point(self, r, x, y, z):
        return self.code.get_gravity_at_point(r, x, y, z)

    def get_potential_at_point(self, r, x, y, z):
        return self.code.get_potential_at_point(r, x, y, z)

    def get_timestep(self):
        return self.code.parameters.timestep

    @property
    def model_time(self):
        return self.code.model_time

    @property
    def particles(self):
        return self.code.particles

    @property
    def total_energy(self):
        return self.code.kinetic_energy + self.code.potential_energy

    @property
    def stop(self):
        return self.code.stop


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
class Gravity(BaseCode):
    def __init__(self, code, particles, eps=0 | units.RSun):
        BaseCode.__init__(self, code, particles, eps)
        self.code.particles.add_particles(self.local_particles)
        self.channel_to_framework = self.code.particles.new_channel_to(
            self.local_particles
        )
        self.channel_from_framework = self.local_particles.new_channel_to(
            self.code.particles
        )
        self.initial_total_energy = self.total_energy


# #BOOKLISTSTOP2# #


# #BOOKLISTSTART3# #
class Hydro(BaseCode):
    def __init__(self, code, particles, eps=0 | units.RSun, dt=None, Rbound=None):
        BaseCode.__init__(self, code, particles, eps)
        self.channel_to_framework = self.code.gas_particles.new_channel_to(
            self.local_particles
        )
        self.channel_from_framework = self.local_particles.new_channel_to(
            self.code.gas_particles
        )
        self.code.gas_particles.add_particles(particles)
        m = self.local_particles.mass.sum()
        l = self.code.gas_particles.position.length()
        if Rbound is None:
            Rbound = 10 * l
        self.code.parameters.periodic_box_size = Rbound
        if dt is None:
            dt = 0.01 * np.sqrt(l**3 / (constants.G * m))
        self.code.parameters.timestep = dt / 8.0
        self.initial_total_energy = self.total_energy

    @property
    def total_energy(self):
        return (
            self.code.kinetic_energy
            + self.code.potential_energy
            + self.code.thermal_energy
        )


# #BOOKLISTSTOP3# #


# #BOOKLISTSTART4# #
def gravity_hydro_bridge(Mprim, Msec, a, ecc, time_end, n_steps, Rgas, Mgas, Ngas):
    stars = new_binary_from_orbital_elements(Mprim, Msec, a, ecc, G=constants.G)
    eps = 1 | units.RSun
    gravity = Gravity(Ph4, stars, eps)

    converter = nbody_system.nbody_to_si(1.0 | units.MSun, Rgas)
    ism = new_plummer_gas_model(Ngas, convert_nbody=converter)
    ism.move_to_center()
    ism = ism.select(lambda r: r.length() < 2 * a, ["position"])
    hydro = Hydro(Fi, ism, eps)
    model_time = 0 | units.Myr
    filename = "gravhydro.hdf5"
    write_set_to_file(stars.savepoint(model_time), filename, "amuse")
    write_set_to_file(ism, filename, "amuse", append_to_file=True)

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,))
    gravhydro.add_system(hydro, (gravity,))
    gravhydro.timestep = 2 * hydro.get_timestep()

    while model_time < time_end:
        orbit = orbital_elements_from_binary(stars, G=constants.G)
        dE_gravity = gravity.initial_total_energy / gravity.total_energy
        dE_hydro = hydro.initial_total_energy / hydro.total_energy
        print(
            (
                "Time:",
                model_time.in_(units.yr),
                "ae=",
                orbit[2].in_(units.au),
                orbit[3],
                "dE=",
                dE_gravity,
                dE_hydro,
            )
        )

        model_time += 10 * gravhydro.timestep
        gravhydro.evolve_model(model_time)
        gravity.copy_to_framework()
        hydro.copy_to_framework()
        write_set_to_file(
            stars.savepoint(model_time), filename, "amuse", append_to_file=True
        )
        write_set_to_file(ism, filename, "amuse", append_to_file=True)
        print("P=", model_time.in_(units.yr), gravity.particles.x.in_(units.au))
    gravity.stop()
    hydro.stop()


# #BOOKLISTSTOP4# #


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-n",
        "--n_steps",
        type=int,
        default=1000,
        help="number of diagnostics time steps",
    )
    parser.add_argument(
        "-N",
        "--Ngas",
        type=int,
        default=1024,
        help="number of gas particles",
    )
    parser.add_argument(
        "--Mprim",
        type=units.MSun,
        dest="Mprim",
        default=3.2 | units.MSun,
        help="Mass of the primary star",
    )
    parser.add_argument(
        "--Msec",
        type=units.MSun,
        dest="Msec",
        default=3.1 | units.MSun,
        help="Mass of the secondary star",
    )
    parser.add_argument(
        "-M",
        type=units.MSun,
        dest="Mgas",
        default=1 | units.MSun,
        help="Mass of the gas",
    )
    parser.add_argument(
        "-R",
        type=units.au,
        dest="Rgas",
        default=10 | units.au,
        help="Size of the gas distribution",
    )
    parser.add_argument(
        "-a",
        type=units.au,
        dest="a",
        default=1 | units.au,
        help="initial orbital separation",
    )
    parser.add_argument(
        "-e",
        dest="ecc",
        type=float,
        default=0.6,
        help="initial orbital eccentricity",
    )
    parser.add_argument(
        "-t",
        type=units.yr,
        dest="time_end",
        default=10000 | units.yr,
        help="end time of the simulation",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    np.random.seed(123)
    gravity_hydro_bridge(**args.__dict__)


if __name__ == "__main__":
    main()
