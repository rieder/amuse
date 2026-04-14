from amuse.units import units
from amuse.community.seba import Seba
from seculartriple_tps import Seculartriple

import numpy as np
import matplotlib.pyplot as plt

# #BOOKLISTSTART1# #
from tres import run_tres, initialize_triple_class
from tres.setup import make_particle_sets


inner_primary_mass = 10.0 | units.MSun
inner_secondary_mass = 8.0 | units.MSun
outer_mass = 5 | units.MSun
inner_semimajor_axis = 1.0 | units.au
outer_semimajor_axis = 12.0 | units.au
inner_eccentricity = 0.5
outer_eccentricity = 0.5
relative_inclination = np.pi / 3
inner_argument_of_pericenter = 0.0
outer_argument_of_pericenter = 0.0
inner_longitude_of_ascending_node = 0.0

stars, bins, correct_params = make_particle_sets(
    inner_primary_mass,
    inner_secondary_mass,
    outer_mass,
    inner_semimajor_axis,
    outer_semimajor_axis,
    inner_eccentricity,
    outer_eccentricity,
    relative_inclination,
    inner_argument_of_pericenter,
    outer_argument_of_pericenter,
    inner_longitude_of_ascending_node,
)

stellar_code = Seba()
secular_code = Seculartriple()


triple = initialize_triple_class(
    stars, bins, correct_params, stellar_code, secular_code
)

triple.stellar_code.parameters.metallicity = 0.02
triple.evolve_model(10 | units.Myr)
triple.print_stellar_system()
# #BOOKLISTSTOP1# #

stellar_code.particles.remove_particles(stars)
stellar_code.stop()
secular_code.stop()
plt.plot(triple.plot_data.times_array.number, triple.plot_data.a_in_array.number)
