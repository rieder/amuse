import logging

import numpy as np

from amuse.units import units, nbody_system
from amuse.datamodel import Particles

# from amuse.community.metisse import Metisse
from amuse_metisse import Metisse
from amuse.community.sse import Sse

from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.plummer import new_plummer_model
from amuse.community.ph4 import Ph4

from amuse.io import write_set_to_file

np.random.seed(127)
logger = logging.getLogger("amuse")
# logger.setLevel(logging.DEBUG)
# logging.basicConfig(level=logging.DEBUG)


def setup_metisse():
    # instance = Metisse(redirection="file")
    instance = Metisse(redirection="none")
    # instance.initialize_code()

    instance.parameters.metallicity_dir = "/Users/rieder/Code/UvA/Toonen/tres3.0/amuse/src/amuse_metisse/data/Hydrogen"
    instance.parameters.metallicity_dir_he = "/Users/rieder/Code/UvA/Toonen/tres3.0/amuse/src/amuse_metisse/data/Helium"

    instance.parameters.wd_mass_scheme = "Modified_mestel"
    instance.parameters.bhns_mass_scheme = "Belczynski2008"
    instance.parameters.initial_metallicity = 0.02

    # instance.commit_parameters()
    return instance


def test_metisse_sun():
    instance = setup_metisse()
    star = Particles(1)
    # star.mass = 4.7893208794726441 | units.MSun
    star.mass = 1.0 | units.MSun

    stars_in_metisse = instance.particles.add_particles(star)
    print(instance.parameters)
    assert stars_in_metisse[0].mass == 1.0 | units.MSun
    # assert stars_in_metisse[0].mass == 4.7893208794726441 | units.MSun
    print(stars_in_metisse[0])
    print("Evolving...")
    instance.evolve_model(1000.0 | units.yr)
    # instance.evolve_one_step(1)
    print(stars_in_metisse[0])
    print("Done")
    instance.stop()


def test_metisse_twostars():
    instance = setup_metisse()
    star = Particles(2)
    star.mass = [0.3, 2.5] | units.MSun

    stars_in_metisse = instance.particles.add_particles(star)
    print(stars_in_metisse[0])
    print("Evolving...")
    instance.evolve_model(1000.0 | units.yr)
    print(stars_in_metisse.stellar_type)
    print("Done")
    instance.stop()

def test_metisse_kroupa():
    instance = setup_metisse()
    number_of_stars = 1000
    star = Particles(number_of_stars)
    star.mass = new_kroupa_mass_distribution(number_of_stars, mass_min=0.3 | units.MSun, mass_max=10.0 | units.MSun)

    stars_in_metisse = instance.particles.add_particles(star)
    print(stars_in_metisse[0])
    print("Evolving...")
    instance.evolve_model(1000.0 | units.yr)
    print(stars_in_metisse)
    print("Done")
    instance.stop()


def evolve_stars_with_metisse(stars, age):
    instance = setup_metisse()
    stars_in_metisse = instance.particles.add_particles(stars)
    instance.evolve_model(age)
    stars_after_evolution = stars_in_metisse.copy()
    instance.stop()
    return stars_after_evolution


def star_cluster_with_metisse(number_of_stars, time_end, time_step, start=0):
    mass = new_kroupa_mass_distribution(
        number_of_stars,
        mass_min=0.75 | units.MSun,
        mass_max=100.0 | units.MSun,
    )
    mass = np.logspace(np.log10(1.0), np.log10(100.0), number_of_stars) | units.MSun
    converter = nbody_system.nbody_to_si(mass.sum(), 3 | units.parsec)
    stars = new_plummer_model(number_of_stars, converter)
    stars.mass = mass

    gravity = Ph4(converter)
    gravity.parameters.epsilon_squared = 0.01 | units.parsec**2
    stars_in_gravity = gravity.particles.add_particles(stars)

    time = 0.0 | units.yr
    i = 0
    time += start * time_step
    i += start
    while time < time_end:
        print(f"Evolving to time: {time}")
        stars_evo = evolve_stars_with_metisse(stars, time)
        # gravity.evolve_model(time)
        # evo_to_gravity = stars_evo.new_channel_to(stars_in_gravity)
        # evo_to_gravity.copy_attributes(["mass"])
        evo_to_model = stars_evo.new_channel_to(stars)
        evo_to_model.copy_attributes(["mass", "luminosity", "stellar_type", "temperature"])
        # grav_to_model = gravity.particles.new_channel_to(stars)
        # grav_to_model.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])
        write_set_to_file(
            stars, f"star_cluster2_metisse_{i:04d}.amuse"
        )
        time += time_step
        i += 1

def evolve_stars_with_sse(stars, age):
    instance = Sse()
    stars_in_sse = instance.particles.add_particles(stars)
    instance.evolve_model(age)
    stars_after_evolution = stars_in_sse.copy()
    instance.stop()
    return stars_after_evolution


def star_cluster_with_metisse_and_sse(number_of_stars, time_end, time_step, start=0):
    mass = np.logspace(np.log10(1.0), np.log10(100.0), number_of_stars) | units.MSun
    stars = Particles(number_of_stars)
    stars.mass = mass

    time = 0.0 | units.yr
    i = 0
    time += start * time_step
    i += start
    while time < time_end:
        print(f"Evolving to time: {time}")
        stars_metisse = evolve_stars_with_metisse(stars, time)
        write_set_to_file(
            stars_metisse, f"stars2_metisse_{i:04d}.amuse"
        )
        # stars_sse = evolve_stars_with_sse(stars, time)
        # write_set_to_file(
        #     stars_sse, f"stars_sse_{i:04d}.amuse"
        # )
        time += time_step
        i += 1

# test_metisse_sun()

# test_metisse_twostars()

# test_metisse_kroupa()
# star_cluster_with_metisse(1000, 100 | units.Myr, 100 | units.kyr, start=0)
star_cluster_with_metisse_and_sse(10000, 100 | units.Myr, 20 | units.kyr, start=0)
