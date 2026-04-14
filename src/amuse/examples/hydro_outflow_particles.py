import os
import argparse
import numpy as np
from amuse.datamodel import Particles, ParticlesSuperset
from amuse.units import units, constants, nbody_system
from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.support.console import set_printing_strategy
from amuse.community.seba import Seba
from amuse.community.fi import Fi
from amuse.io import write_set_to_file


set_printing_strategy(
    "custom",
    preferred_units=[units.MSun, units.au, units.Myr],
    precision=5,
    prefix="",
    separator=" [",
    suffix="]",
)


# #BOOKLISTSTART1# #
def new_sph_particles_from_stellar_wind(stars, mass_per_gas_particle):
    new_sph = Particles(0)
    for i, star in enumerate(stars):
        number_of_gas = int(-star.Mwind / mass_per_gas_particle)
        if number_of_gas == 0:
            continue
        total_mass_of_gas = mass_per_gas_particle * number_of_gas
        star.Mwind += total_mass_of_gas
        new_gas = Particles(number_of_gas)
        new_gas.mass = mass_per_gas_particle
        new_gas.h_smooth = 0.0 | units.parsec

        dx, dy, dz = uniform_unit_sphere(number_of_gas).make_xyz()
        new_gas.x = star.x + (dx * star.radius)
        new_gas.y = star.y + (dy * star.radius)
        new_gas.z = star.z + (dz * star.radius)
        for j, gas_particle in enumerate(new_gas):
            r = gas_particle.position - star.position
            r = r / r.length()
            v_wind = (
                constants.G
                * star.mass
                / (gas_particle.position - star.position).length()
            ).sqrt()
            gas_particle.u = 0.5 * (v_wind) ** 2
            gas_particle.vx = star.vx + r[0] * star.terminal_wind_velocity
            gas_particle.vy = star.vy + r[1] * star.terminal_wind_velocity
            gas_particle.vz = star.vz + r[2] * star.terminal_wind_velocity
        new_sph.add_particles(new_gas)
    return new_sph


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
def v_terminal_teff(star):
    t4 = (np.log10(star.temperature.value_in(units.K)) - 4.0).clip(0.0, 1.0)
    return (30 | units.km / units.s) + ((4000 | units.km / units.s) * t4)


# #BOOKLISTSTOP2# #


def hydro_outflow_particles():
    stars = Particles(2)
    stars.mass = (9.5, 10) | units.MSun
    stars[0].position = (1, 0, 0) | units.au
    stars[0].velocity = (0, 0, 0) | units.kms
    stars[1].position = (0, 0, 0) | units.au
    stars[1].velocity = (0, 0, 0) | units.kms
    stars.move_to_center()

    a = stars.position.length().amax()
    vc = constants.G * stars.mass.sum() / a
    stellar = Seba()
    stellar.particles.add_particles(stars)
    stellar_to_framework = stellar.particles.new_channel_to(stars)
    stellar.evolve_model(26 | units.Myr)
    stellar_to_framework.copy_attributes(["mass", "radius", "temperature"])
    dt = 0.1 | units.Myr
    stellar.evolve_model((26 | units.Myr) + dt)
    stars.dmdt = (stellar.particles.mass - stars.mass) / dt
    stars.Mwind = 0 | units.MSun
    stars.terminal_wind_velocity = v_terminal_teff(stars)
    stellar.stop()
    dt = 0.1 | units.day
    mgas = 0.1 * abs(stars.dmdt.sum() * dt)

    converter = nbody_system.nbody_to_si(1 | units.MSun, a)
    bodies = Particles(0)
    bodies.mass = mgas
    bodies.position = (0, 0, 0) | units.au
    bodies.velocity = (0, 0, 0) | units.kms
    bodies.u = 0 | units.m**2 * units.s**-2
    bodies.h_smooth = 0.01 * a

    hydro = Fi(converter, redirection="none")
    if len(bodies) > 0:
        hydro.gas_particles.add_particles(bodies)
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.timestep = dt
    hydro.parameters.periodic_box_size = 1000 * a
    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)

    moving_bodies = ParticlesSuperset([stars, bodies])
    filename = "hydro_outflow_particles.amuse"
    istep = 0
    while hydro.model_time < 10 | units.yr:
        stars.Mwind += stars.dmdt * dt
        new_sph = new_sph_particles_from_stellar_wind(stars, mgas)
        if len(new_sph) > 0:
            bodies.add_particles(new_sph)
            bodies.synchronize_to(hydro.gas_particles)
        print(f"time={hydro.model_time} Ngas={len(bodies)} {mgas * len(bodies)}")
        if len(bodies) > 100:
            hydro.evolve_model(hydro.model_time + dt)
            hydro_to_framework.copy()
            if istep % 10 == 0:
                write_set_to_file(moving_bodies, filename, append_to_file=True)
            istep += 1
    hydro.stop()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    return parser


def main(**kwargs):
    hydro_outflow_particles(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
