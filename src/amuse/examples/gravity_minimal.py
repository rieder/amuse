"""
Minimal example for gravity calculations in AMUSE
"""

from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.hermite import Hermite


def gravity_minimal(number_of_stars=100, time_end=1 | nbody_system.time):
    """
    Initialises a Plummer sphere of length 'number_of_stars' and integrates until
    'time_end' with the Hermite integrator
    """
    stars = new_plummer_model(number_of_stars)

    gravity = Hermite()
    gravity.particles.add_particles(stars)

    energy_total_initial = gravity.kinetic_energy + gravity.potential_energy

    gravity.evolve_model(time_end)

    energy_total = gravity.kinetic_energy + gravity.potential_energy
    energy_error = (energy_total_initial - energy_total) / energy_total_initial

    print(
        f"model time = {gravity.model_time} "
        f"total mass = {stars.total_mass()} "
        f"total energy = {energy_total} "
        f"relative energy error = {energy_error}"
    )

    gravity.stop()


def main(**kwargs):
    gravity_minimal(**kwargs)


if __name__ == "__main__":
    main(number_of_stars=100, time_end=1 | nbody_system.time)
