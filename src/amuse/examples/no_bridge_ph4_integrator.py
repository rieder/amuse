from amuse.units import nbody_system, units
from amuse.units.quantities import zero
from amuse.io import write_set_to_file
from amuse.couple import bridge
from amuse.community.ph4 import Ph4
import numpy as np

from new_system_of_sun_earth_and_moon import new_system_of_sun_earth_and_moon

def main():
    filename = "../data/SunAndEarthAndMoon_ISS.amuse"
    ss = new_system_of_sun_earth_and_moon()
    
    converter = nbody_system.nbody_to_si(ss.mass.sum(), ss.position.length())
    gravity = Ph4(converter)
    gravity.parameters.timestep_parameter = 0.15
    gravity.particles.add_particle(ss)

    channel_to_framework = gravity.particles.new_channel_to(ss)

    time = 0 | units.yr
    write_set_to_file(ss, filename,
                      timestamp=time, overwrite_file=True)

    energy_total_init = gravity.kinetic_energy + gravity.potential_energy
    energy_total_prev = energy_total_init

    gravity.timestep = 1 | units.day
    time_step = 0.1 | units.yr
    time_end = 20 | units.yr
    ddE_max = 0
    while time < time_end:
        time += time_step
        gravity.evolve_model(time)

        # Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_to_framework.copy()
        write_set_to_file(ss, filename, 
                          timestamp=time, append_to_file=True, overwrite_file=False)

        energy_kinetic = gravity.kinetic_energy
        energy_potential = gravity.potential_energy
        energy_total = energy_kinetic + energy_potential
        ddE_max = max(ddE_max, np.abs((energy_total_prev-energy_total)/energy_total))
        print(
            f"T={time} "
            f"E= {energy_total} Q= {energy_kinetic/energy_potential} "
            f"dE= {(energy_total_init-energy_total)/energy_total} "
            f"ddE= {(energy_total_prev-energy_total)/energy_total} "
            f"ddE_Max= {ddE_max}"
        )
        energy_total_prev = energy_total
    gravity.stop()


if __name__ in ('__main__', '__plot__'):
    main()
