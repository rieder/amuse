"""
Example to merge two stars using SPH
"""
import argparse
import numpy as np
import matplotlib.pyplot as plt

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particle, ParticlesSuperset
from amuse.community.evtwin import Evtwin
from amuse.community.gadget2 import Gadget2
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH


def return_evolved_star_hydro(mass, time, Nsph):
    star = Particle(mass=mass)
    stellar = Evtwin()
    star = stellar.particles.add_particle(star)
    stellar.evolve_model(time)
    Nsph = Nsph * int(mass.value_in(units.MSun))
    star_in_sph = convert_stellar_model_to_SPH(star, Nsph).gas_particles
    stellar.stop()
    return star_in_sph


def merge_two_stars_sph(
    mass_primary,
    mass_secondary,
    time_collision,
    Nsph,
    opening_angle,
):
    primary_in_sph = return_evolved_star_hydro(
        mass_primary, time_collision, int(Nsph*mass_primary/(1. | units.MSun))
    )
    # primary_in_sph = relax_sph_realization(primary_in_sph)
    secondary_in_sph = return_evolved_star_hydro(
        mass_secondary,
        time_collision,
        int(Nsph*mass_secondary/(1. | units.MSun))
    )
    # secondary_in_sph = relax_sph_realization(secondary_in_sph)
    radius = primary_in_sph.x.max() + secondary_in_sph.x.max()
    total_mass = primary_in_sph.mass.sum() + secondary_in_sph.mass.sum()
    secondary_in_sph.x += 0.8*radius
    secondary_in_sph.y += 0.6*radius
    secondary_in_sph.vx -= (constants.G*total_mass/radius).sqrt()

    converter = nbody_system.nbody_to_si(mass_primary, 1.0 | units.au)
    hydro = Gadget2(converter, number_of_workers=4)
    hydro.parameters.opening_angle = opening_angle

    # print("Opening criterion:", hydro.parameters.opening_angle)
    # print("Opening criterion:", hydro.get_gdgop())
    # print("Opening criterion:", hydro.parameters.opening_angle)
    # print("Opening criterion:", hydro.get_gdgop())

    hydro.gas_particles.add_particles(primary_in_sph)
    hydro.gas_particles.add_particles(secondary_in_sph)
    hydro.evolve_model(2.0 | units.hour)
    hydro.gas_particles.new_channel_to(primary_in_sph).copy()
    hydro.gas_particles.new_channel_to(secondary_in_sph).copy()
    hydro.stop()
    return primary_in_sph, secondary_in_sph


def relax_sph_realization(sph_star):
    dynamical_timescale = sph_star.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1 | units.RSun)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.gas_particles.add_particles(sph_star)

    # to_hydro = sph_star.new_channel_to(hydro.gas_particles)
    to_framework = hydro.gas_particles.new_channel_to(sph_star)

    ts_factor = 2.5
    t_end = ts_factor * sph_star.dynamical_timescale(mass_fraction=0.9)
    n_steps = ts_factor * 100
    velocity_damp_factor = 1.0 - (ts_factor*2*np.pi)/n_steps
    dt = t_end/float(n_steps)
    time = 0 | units.day
    while time < t_end:
        time += dt
        hydro.evolve_model(time)
        hydro.gas_particles.velocity = (
            velocity_damp_factor
            * hydro.gas_particles.velocity
        )
    to_framework.copy()
    hydro.stop()
    return sph_star


def merge_and_plot_distribution(
    mass_primary,
    mass_secondary,
    time_collision,
    Nsph,
    opening_angle,
    c=None,
    label=None
):
    p, s = merge_two_stars_sph(
        mass_primary,
        mass_secondary,
        time_collision,
        Nsph,
        opening_angle,
    )
    merger = ParticlesSuperset([p, s])
    com = merger.center_of_mass()
    merger.r = (
        (merger.x-com[0])**2
        + (merger.y-com[1])**2
        + (merger.z-com[2])**2
    ).sqrt()
    merger = merger.sorted_by_attributes("r")
    n = []
    m = 0
    mi = 1.0*(mass_primary+mass_secondary).value_in(units.MSun)/len(merger)
    for i in range(len(merger.r)):
        m += mi
        n.append(m)
    plt.plot(merger.r.value_in(units.RSun), n, c=c, label=label)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-M", "--mass_primary", type=units.MSun, default=10 | units.MSun,
        help="Mass of the primary star"
    )
    parser.add_argument(
        "-m", "--mass_secondary", type=units.MSun, default=1 | units.MSun,
        help="Mass of the secondary star"
    )
    # specifying more numbers will result in more simulations
    parser.add_argument(
        "-N", "--number_of_sph_particles", type=int, nargs="+",
        default=[10, 100, 1000, 10000],
        help="number of sph particles"
    )
    parser.add_argument(
        "-t",
        "--time_collision",
        type=units.Myr,
        default=0.01 | units.Myr,
        help="end time of the simulation"
    )
    return parser


def merge_two_stars_sph_convergence(
    mass_primary=10 | units.MSun,
    mass_secondary=1 | units.MSun,
    time_collision=0.01 | units.Myr,
    number_of_sph_particles=None,
):
    if number_of_sph_particles is None:
        number_of_sph_particles = [10, 100, 1000, 10000]
    plt.figure(figsize=(8, 6))
    plt.xlabel(r"R [R$_\odot$]")
    plt.ylabel("$M_{<r}$")
    ax = plt.gca()
    ax.minorticks_on()  # switch on the minor ticks
    ax.locator_params(nbins=3)
    for n in number_of_sph_particles:
        merge_and_plot_distribution(
            mass_primary, mass_secondary, time_collision, n, 0.5,
            label=r"$N_{sph} = " + f"{n}" + r"/M_\odot$"
        )
    plt.xlim(0, 20)
    plt.ylim(0, 12)
    plt.legend(loc=4, fontsize=24)
    plt.savefig("stellar_merger_convergence.pdf")


def main():
    args = new_argument_parser().parse_args()
    merge_two_stars_sph_convergence(**args.__dict__)


if __name__ == "__main__":
    main()
