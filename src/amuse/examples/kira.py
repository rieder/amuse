import logging
import argparse
import random

import numpy as np
from amuse.community.kepler import Kepler
from amuse.community.smalln import Smalln
from amuse.community.hermite import Hermite
from amuse.community.seba import Seba
from amuse.datamodel import Particles
from amuse.units.quantities import as_vector_quantity
from amuse.couple import encounters
from amuse.units import constants, units, quantities, nbody_system
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.plummer import new_plummer_model
from amuse.support.console import set_printing_strategy

from matplotlib import pyplot as plt


def new_smalln(converter):
    result = Smalln(converter)
    result.parameters.timestep_parameter = 0.1
    result.parameters.cm_index = 2001
    return result


def new_kepler(converter):
    kepler = Kepler(converter)
    kepler.initialize_code()
    kepler.set_longitudinal_unit_vector(1.0, 0.0, 0.0)
    kepler.set_transverse_unit_vector(0.0, 1.0, 0)
    return kepler


def new_binary_orbit(mass1, mass2, semi_major_axis, eccentricity=0, keyoffset=1):
    total_mass = mass1 + mass2
    mass_fraction_particle_1 = mass1 / (total_mass)

    binary = Particles(2)
    binary[0].mass = mass1
    binary[1].mass = mass2

    mu = constants.G * total_mass

    velocity_perihelion = np.sqrt(
        mu / semi_major_axis * ((1.0 + eccentricity) / (1.0 - eccentricity))
    )
    radius_perihelion = semi_major_axis * (1.0 - eccentricity)
    print(velocity_perihelion)

    binary[0].position = (
        (1.0 - mass_fraction_particle_1) * radius_perihelion * [1.0, 0.0, 0.0]
    )
    binary[1].position = -(
        mass_fraction_particle_1 * radius_perihelion * [1.0, 0.0, 0.0]
    )

    binary[0].velocity = (
        (1.0 - mass_fraction_particle_1) * velocity_perihelion * [0.0, 1.0, 0.0]
    )
    binary[1].velocity = -(
        mass_fraction_particle_1 * velocity_perihelion * [0.0, 1.0, 0.0]
    )

    return binary


# see Eggleton 2006 Equation 1.6.3 (2006epbm.book.....E)
def random_semimajor_axis_PPE(
    Mprim, Msec, P_min=10.0 | units.day, P_max=100.0 | units.yr
):

    Pmax = P_max.value_in(units.day)
    Pmin = P_min.value_in(units.day)
    mpf = (Mprim.value_in(units.MSun) ** 2.5) / 5.0e4
    rnd_max = (Pmax * mpf) ** (1.0 / 3.3) / (1 + (Pmin * mpf) ** (1.0 / 3.3))
    rnd_min = (Pmin * mpf) ** (1.0 / 3.3) / (1 + (Pmax * mpf) ** (1.0 / 3.3))
    rnd_max = min(rnd_max, 1)
    rnd = np.random.uniform(rnd_min, rnd_max, 1)
    Porb = ((rnd / (1.0 - rnd)) ** 3.3) / mpf | units.day
    Mtot = Mprim + Msec
    a = ((constants.G * Mtot) * (Porb / (2 * np.pi)) ** 2) ** (1.0 / 3.0)
    return a


def make_secondaries(center_of_masses, Nbin):

    resulting_binaries = Particles()
    singles_in_binaries = Particles()
    binaries = center_of_masses.random_sample(Nbin)
    mmin = center_of_masses.mass.min()
    for bi in binaries:
        mp = bi.mass
        ms = (
            np.random.uniform(mmin.value_in(units.MSun), mp.value_in(units.MSun))
            | units.MSun
        )
        a = random_semimajor_axis_PPE(mp, ms)
        e = np.sqrt(np.random.random())

        nb = new_binary_orbit(mp, ms, a, e)
        nb.position += bi.position
        nb.velocity += bi.velocity
        nb = singles_in_binaries.add_particles(nb)
        nb.radius = 0.01 * a

        bi.radius = 3 * a
        binary_particle = bi.copy()
        binary_particle.child1 = nb[0]
        binary_particle.child2 = nb[1]
        binary_particle.semi_major_axis = a
        binary_particle.eccentricity = e
        resulting_binaries.add_particle(binary_particle)

    single_stars = center_of_masses - binaries
    return single_stars, resulting_binaries, singles_in_binaries


def calculate_orbital_elementss(bi, converter):
    kep = new_kepler(converter)
    comp1 = bi.child1
    comp2 = bi.child2
    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a, e = kep.get_elements()
    kep.stop()
    return a, e


# #BOOKLISTSTART# #


def resolve_changed_binaries(stopping_condition, stellar, converter):
    new_binaries = stopping_condition.particles(0)
    for bi in new_binaries:
        print(("add new binary:", bi))
        a, e = calculate_orbital_elementss(bi, converter)
        bi.semi_major_axis = a
        bi.eccentricity = e
        stellar.binaries.add_particle(bi)
        print(("new binary parameters", a, e))
        print(bi)

    lost_binaries = stopping_condition.particles(1)
    for bi in lost_binaries:
        print(("remove old binary:", bi.key))
        stellar.binaries.remove_particle(bi)

    changed_binaries = stopping_condition.particles(2)
    for bi in changed_binaries:
        bs = bi.as_particle_in_set(stellar.binaries)
        a, e = calculate_orbital_elementss(bi, converter)
        bs.semi_major_axis = a
        bs.eccentricity = e
        print(("Modified binary parameters", a, e))
        print(bs)


# #BOOKLISTSTOP# #


def update_dynamical_binaries_from_stellar(stellar, multiples_code, converter):
    kep = new_kepler(converter)

    # THIS NEEDS TO BE CHECKED!
    print("++++++++++++ THIS NEEDS TO BE CHECKED ++++++++++++++++++++")

    print(("Number of binaries=", len(stellar.binaries)))
    for bi in stellar.binaries:
        bs = bi.as_particle_in_set(multiples_code.binaries)
        total_mass = bi.child1.mass + bi.child2.mass
        kep.initialize_from_elements(total_mass, bi.semi_major_axis, bi.eccentricity)
        rel_position = as_vector_quantity(kep.get_separation_vector())
        rel_velocity = as_vector_quantity(kep.get_velocity_vector())
        mu = bi.child1.mass / total_mass
        bs.child1.position = mu * rel_position
        bs.child2.position = -(1 - mu) * rel_position
        bs.child1.velocity = mu * rel_velocity
        bs.child2.velocity = -(1 - mu) * rel_velocity
        print(
            (
                "semi_major_axis=",
                bi.semi_major_axis,
                total_mass,
                bi.child1.mass,
                bi.child2.mass,
                bi.eccentricity,
            )
        )
    kep.stop()


def kira(
    time_end=10.0 | units.Myr,
    number_of_particles=100,
    radius=1.0 | units.parsec,
    number_of_binaries=50,
    seed=2,
):
    if seed >= 0:
        np.random.seed(seed)
        # This is only for random.sample, which apparently does not use numpy
        random.seed(seed)
    logging.basicConfig(level=logging.ERROR)

    mass = new_salpeter_mass_distribution(number_of_particles, mass_min=10 | units.MSun)
    converter = nbody_system.nbody_to_si(mass.sum(), radius)
    code = Hermite(converter)
    stars = new_plummer_model(number_of_particles, convert_nbody=converter)
    stars.mass = mass
    stars.radius = 0.01 / len(stars) | radius.unit

    single_stars, binary_stars, singles_in_binaries = make_secondaries(
        stars, number_of_binaries
    )
    print(binary_stars)

    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    figure = plt.figure(figsize=(8, 6))
    ax = figure.add_subplot(1, 1, 1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.xaxis._autolabelpos = True
    ax.yaxis._autolabelpos = True
    ax.set_xscale("log")
    ax.set_xlim(1.0e-3, 1.0e3)
    ax.set_ylim(-0.04, 1.04)

    # Just to make the point in the legend at the top.
    ax.scatter(
        [1], [10], c=color[0], s=100, lw=0, alpha=1.0, label="Initial binary parameters"
    )

    a_init = binary_stars.semi_major_axis
    e_init = binary_stars.eccentricity
    # ax.scatter(binary_stars.semi_major_axis.value_in(units.au),
    #            binary_stars.eccentricity, c=color[0], s=50, lw=0,
    #            alpha=0.5,
    #            label="Initial binary parameters")

    stellar = Seba()
    stellar.particles.add_particles(single_stars)
    stellar.particles.add_particles(singles_in_binaries)
    stellar.binaries.add_particles(binary_stars)
    channel_to_stars = stellar.particles.new_channel_to(stars)

    encounter_code = encounters.HandleEncounter(
        kepler_code=new_kepler(converter),
        resolve_collision_code=new_smalln(converter),
        interaction_over_code=None,
        G=constants.G,
    )
    multiples_code = encounters.Multiples(
        gravity_code=code, handle_encounter_code=encounter_code, G=constants.G
    )
    multiples_code.particles.add_particles((stars - binary_stars).copy())
    multiples_code.singles_in_binaries.add_particles(singles_in_binaries)
    multiples_code.binaries.add_particles(binary_stars)
    multiples_code.commit_particles()
    channel_from_stars_to_particles = stellar.particles.new_channel_to(
        multiples_code.particles
    )

    stopping_condition = multiples_code.stopping_conditions.binaries_change_detection
    stopping_condition.enable()

    ax.scatter(
        stellar.binaries.semi_major_axis.value_in(units.au),
        stellar.binaries.eccentricity,
        c=color[0],
        s=200,
        lw=0,
        alpha=0.5,
        label="After binary evolution",
    )

    t = quantities.linspace(0 * time_end, time_end, 11)
    for ti in t:
        print(
            f"t, Energy= {ti} {multiples_code.particles.mass.sum()} "
            f"{multiples_code.get_total_energy()}"
        )
        multiples_code.evolve_model(ti)
        print(
            f"at t={multiples_code.model_time} "
            f"Nmultiples: {len(multiples_code.multiples)}"
        )

        if stopping_condition.is_set():
            resolve_changed_binaries(stopping_condition, stellar, converter)

        stellar.evolve_model(ti)
        channel_from_stars_to_particles.copy_attributes(["mass", "radius"])
        update_dynamical_binaries_from_stellar(stellar, multiples_code, converter)

        print(
            f"Lagrangian radii: {multiples_code.all_singles.LagrangianRadii(converter)}"
        )
        print(f"MC.particles {multiples_code.particles}")
        print(
            f"Lagrangian radii: {multiples_code.particles.LagrangianRadii(converter)}"
        )
        print(f"t, Energy= {ti, multiples_code.get_total_energy()}")

    ax.scatter(a_init.value_in(units.au), e_init, c=color[0], s=100, lw=0, alpha=1.0)

    print(f"Stellar type: {stellar.particles.stellar_type}")
    ax.scatter(
        stellar.binaries.semi_major_axis.value_in(units.au),
        stellar.binaries.eccentricity,
        c=color[1],
        lw=0,
        s=50,
        label="Binary evolution with dynamics",
    )

    ax.set_xlabel(r"$\log_{10}(a/R_\odot)$")
    ax.set_ylabel("e")
    ax.legend(loc="lower right")

    save_file = "kira_a_vs_e.pdf"
    plt.savefig(save_file)
    print(f"\nSaved figure in file {save_file} \n")

    stellar.stop()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-t", "--time_end", type=units.Myr, default=10.0 | units.Myr)
    parser.add_argument("-R", "--radius", type=units.parsec, default=1.0 | units.parsec)
    parser.add_argument("-N", "--number_of_particles", type=int, default=100)
    parser.add_argument("--Nbin", "--number_of_binaries", type=int, default=50)
    parser.add_argument("--seed", type=int, default=2)
    return parser


def main(**kwargs):
    kira(**kwargs)


if __name__ == "__main__":
    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.parsec, units.Myr],
        precision=4,
        prefix="",
        separator=" [",
        suffix="]",
    )

    arguments = new_argument_parser().parse_args()
    kira(**arguments.__dict__)
