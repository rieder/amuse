import logging

import numpy as np

from amuse.units import units
from amuse.datamodel import Particles

# from amuse.community.metisse import Metisse
from amuse_metisse import Metisse

logger = logging.getLogger("amuse")
# logger.setLevel(logging.DEBUG)
# logging.basicConfig(level=logging.DEBUG)


def setup_metisse():
    # instance = Metisse(redirection="file")
    instance = Metisse(redirection="none")
    # instance.initialize_code()

    instance.parameters.metallicity_dir = "/Users/rieder/Code/UvA/Toonen/tres3.0/amuse/src/amuse_metisse/data/hydrogen"
    instance.parameters.metallicity_dir_he = "/Users/rieder/Code/UvA/Toonen/tres3.0/amuse/src/amuse_metisse/data/helium"

    instance.parameters.wd_mass_scheme = "Modified_mestel"
    instance.parameters.bhns_mass_scheme = "Belczynski2008"
    instance.parameters.initial_metallicity = 0.02

    # instance.commit_parameters()
    return instance


def test_metisse_sun():
    instance = setup_metisse()
    star = Particles(1)
    star.mass = 4.7893208794726441 | units.MSun
    # star.mass = 1.0 | units.MSun

    stars_in_metisse = instance.particles.add_particles(star)
    print(instance.parameters)
    # assert stars_in_metisse[0].mass == 1.0 | units.MSun
    assert stars_in_metisse[0].mass == 4.7893208794726441 | units.MSun
    print(stars_in_metisse[0])
    print("Evolving...")
    # instance.evolve_model(1000.0 | units.yr)
    instance.evolve_one_step(1)
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
    from amuse.ic.kroupa import new_kroupa_mass_distribution
    np.random.seed(127)
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


test_metisse_sun()

# test_metisse_twostars()

# test_metisse_kroupa()
