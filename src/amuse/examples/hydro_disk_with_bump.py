"""
Simulate the radiative and hydrodynamial evolution of a disk with a bump
around a single star
"""

import argparse
from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles
from amuse.io import write_set_to_file
from amuse.community.simplex.interface import SimpleXSplitSet
from amuse.community.gadget2 import Gadget2
from amuse.community.seba import Seba
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.support.console import set_printing_strategy

set_printing_strategy(
    "custom",
    preferred_units=[units.MSun, units.au, units.Myr],
    precision=12,
    prefix="",
    separator=" [",
    suffix="]",
)


def mu(X=None, Y=0.25, Z=0.02, x_ion=0.1):
    """
    Compute the mean molecular weight in kg (the average weight of
    particles in a gas) X, Y, and Z are the mass fractions of
    Hydrogen, of Helium, and of metals, respectively.  x_ion is the
    ionisation fraction (0 < x_ion < 1), 1 means fully ionised.
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception(
            "Error in calculating mu: mass " + "fractions do not sum to 1.0"
        )
    return constants.proton_mass / (
        X * (1.0 + x_ion) + Y * (1.0 + 2.0 * x_ion) / 4.0 + Z * x_ion / 2.0
    )


class RadHydro:
    def __init__(self, rad, hydro, star, disk):

        self.time = 0 | units.day
        self.star = star
        self.disk = disk

        disk.r2 = disk.x**2 + disk.y**2
        disk = disk.sorted_by_attributes("r2")
        radius_max = disk.r2.max().sqrt()
        print("MaxR=", radius_max.in_(units.au))

        self.hydro = hydro(
            nbody_system.nbody_to_si(self.disk.mass.sum(), radius_max),
            number_of_workers=4,
        )
        self.hydro.parameters.epsilon_squared = (10 | units.au) ** 2

        self.hydro.gas_particles.add_particles(self.disk)
        self.hydro.gas_particles.new_channel_to(self.disk)
        self.hydro.dm_particles.add_particles(self.star)
        self.hydro.dm_particles.new_channel_to(self.star)

        self.hydro_to_star = self.hydro.dm_particles.new_channel_to(self.star)
        self.hydro_to_disk = self.hydro.gas_particles.new_channel_to(self.disk)
        self.star_to_hydro = self.star.new_channel_to(self.hydro.dm_particles)
        self.disk_to_hydro = self.disk.new_channel_to(self.hydro.gas_particles)

        self.hydro.evolve_model(1 | units.s)
        self.hydro_to_star.copy()
        self.hydro_to_disk.copy()

        self.rad = rad()
        for si in self.star:
            if si.mass >= 5 | units.MSun:
                self.rad.src_particles.add_particle(si)
        self.rad.gas_particles.add_particles(self.disk)
        self.rad.parameters.box_size = 2.01 * radius_max
        self.rad.parameters.timestep = 1 | units.day
        self.rad.set_source_Teff(star.temperature)

        self.rad_to_disk = self.rad.gas_particles.new_channel_to(
            self.disk,
            attributes=["xion", "u"],
        )
        self.star_to_rad = self.star.new_channel_to(
            self.rad.src_particles,
            attributes=["x", "y", "z"],
        )
        self.disk_to_rad = self.disk.new_channel_to(
            self.rad.gas_particles,
            attributes=["x", "y", "z"],
        )
        # self.rad_to_disk.copy()
        # self.star_to_rad.copy()
        # self.disk_to_rad.copy()

        self.rad.stop()
        self.index = 0

    def write_file(self):
        self.index += 1
        filename = f"hydro_disk_with_bump_i{self.index:04}.amuse"
        write_set_to_file(self.star, filename, overwrite_file=True)
        write_set_to_file(self.disk, filename, append_to_file=True)

    def evolve_model(self, model_time):
        dt = model_time - self.time
        self.old_time = self.time
        self.time += dt / 2.0

        # self.disk_to_rad.copy()
        # self.star_to_rad.copy()
        # self.rad.evolve_model(self.time)
        # self.rad_to_disk.copy()

        self.time += dt / 2.0

        self.disk_to_hydro.copy()
        self.star_to_hydro.copy()
        self.hydro.evolve_model(self.time)
        self.hydro_to_disk.copy()
        self.hydro_to_star.copy()

        print("RT done at time:", self.time.in_(units.day))

    def print_diagnostics(self):
        umin = self.disk.u.min()
        umean = self.disk.u.mean()
        umax = self.disk.u.max()
        Tmin = mu() / constants.kB * umax
        Tmean = mu() / constants.kB * umean
        Tmax = mu() / constants.kB * umin

        print("Time=", self.time.in_(units.day))
        print(
            f"Ionization: {self.disk.xion.min()} {self.disk.xion.mean()}"
            f"self.disk.xion.max()"
        )
        print("Intenal energy:", umin, umean, umax)
        print("Temperature:", Tmin, Tmean, Tmax)
        print(
            "Density:",
            self.disk.density.min().in_(units.amu / units.cm**3),
            self.disk.density.mean().in_(units.amu / units.cm**3),
            self.disk.density.max().in_(units.amu / units.cm**3),
        )
        print("scaleheight:", abs(self.disk.z.value_in(units.au)).mean())

    def stop(self):
        self.hydro.stop()


# #BOOKLISTSTART1# #
def new_disk_with_bump(
    mass_star=10 | units.MSun,
    n_disk=100,
    mass_disk=1.0 | units.MSun,
    radius_min=1.0 | units.au,
    radius_max=100.0 | units.au,
    mass_bump=0.1 | units.MSun,
    radius_bump=5.0 | units.au,
    dist_bump=10 | units.au,
):

    converter = nbody_system.nbody_to_si(mass_star, radius_min)
    disk = ProtoPlanetaryDisk(
        n_disk,
        convert_nbody=converter,
        densitypower=1.5,
        Rmin=1,
        Rmax=radius_max / radius_min,
        q_out=1.0,
        discfraction=mass_disk / mass_star,
    ).result
    com = disk.center_of_mass()

    # determine bump's local velocity

    inner_particles = disk.select(
        lambda r: (com - r).length() < dist_bump, ["position"]
    )
    M_inner = mass_star + inner_particles.mass.sum()
    v_circ = (
        (constants.G * M_inner * (2.0 / dist_bump - 1.0 / dist_bump)).sqrt().value_in(units.kms)
    )

    # initialize bump

    Nbump = int(n_disk * mass_bump / mass_disk)
    bump = new_plummer_gas_model(
        Nbump, convert_nbody=nbody_system.nbody_to_si(mass_bump, radius_bump)
    )
    bump.x += dist_bump
    bump.velocity += [0, v_circ, 0] | units.kms

    disk.add_particles(bump)
    disk.move_to_center()
    return disk


# #BOOKLISTSTOP1# #


def evolve_star(mass_star, tstar):
    stars = Particles(1)
    stars = Particles(2)
    stars[0].mass = mass_star
    stars[1].mass = 0.1 * mass_star
    stellar = Seba()
    stellar.particles.add_particle(stars)
    stellar.evolve_model(tstar)
    stars.mass = stellar.particles.mass
    stars.position = (0, 0, 0) | units.au
    stars.velocity = (0, 0, 0) | units.kms
    stars.luminosity = stellar.particles.luminosity / (20.0 | units.eV)
    stars.temperature = stellar.particles.temperature
    stars.flux = stars.luminosity
    stars.rho = 1.0 | (units.g / units.cm**3)
    stars.xion = 0.0  # ionization_fraction
    stars.u = (9.0 | units.kms) ** 2  # internal_energy
    print(stars)

    if len(stars) > 1:
        stars[1].x = 50 | units.au
        vc = 1.0 * (constants.G * stars.mass.sum() / (100.0 | units.au)).sqrt()
        stars[1].vy += vc

    stellar.stop()
    return stars


def hydro_disk_with_bump(
    mass_star=10 | units.MSun,
    n_disk=100,
    mass_disk=1.0 | units.MSun,
    radius_min=1.0 | units.au,
    radius_max=100.0 | units.au,
    mass_bump=0.1 | units.MSun,
    radius_bump=5.0 | units.au,
    dist_bump=10 | units.au,
    time_end=10 | units.yr,
    n_steps=10,
):

    star = evolve_star(mass_star, time_end)

    disk = new_disk_with_bump(
        star[0].mass,
        n_disk,
        mass_disk,
        radius_min,
        radius_max,
        mass_bump,
        radius_bump,
        dist_bump,
    )

    radhydro = RadHydro(SimpleXSplitSet, Gadget2, star, disk)
    radhydro.write_file()

    time_step = time_end / float(n_steps)
    print(f"time_step={time_step.in_(units.day)}")
    time_model = 0 | units.day
    while time_model < time_end:
        time_model += time_step
        radhydro.evolve_model(time_model)
        radhydro.print_diagnostics()
        radhydro.write_file()
    radhydro.stop()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-N",
        "--n_disk",
        type=int,
        default=10000,
        help="number of particles in disk",
    )
    parser.add_argument(
        "-t",
        "--time_end",
        type=units.yr,
        default=2000.0 | units.yr,
        help="radiation time",
    )
    parser.add_argument(
        "-n", "--n_steps", type=int, default=100, help="number of steps"
    )
    parser.add_argument(
        "--mass_star",
        type=units.MSun,
        default=10 | units.MSun,
        help="Mass of the central star",
    )
    parser.add_argument(
        "--mass_disk",
        type=units.MSun,
        default=1 | units.MSun,
        help="Mass of the disk",
    )
    parser.add_argument(
        "-r",
        "--radius_min",
        type=units.au,
        default=10 | units.au,
        help="inner disk radius",
    )
    parser.add_argument(
        "-R",
        "--radius_max",
        type=units.au,
        default=100 | units.au,
        help="outer disk radius",
    )
    parser.add_argument(
        "--mass_bump",
        type=units.MSun,
        default=0.5 | units.MSun,
        help="bump mass",
    )
    parser.add_argument(
        "--radius_bump",
        type=units.au,
        default=5 | units.au,
        help="bump radius",
    )
    parser.add_argument(
        "-a",
        "--dist_bump",
        type=units.au,
        default=50 | units.au,
        help="distance of bump from star",
    )
    return parser


def main():
    arguments = new_argument_parser().parse_args()
    hydro_disk_with_bump(**arguments.__dict__)


if __name__ == "__main__":
    main()
