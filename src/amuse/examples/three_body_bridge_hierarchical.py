import numpy as np
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.units.quantities import zero
from amuse.community.ph4 import Ph4
from amuse.io import write_set_to_file
from amuse.couple import bridge

from new_system_of_sun_earth_and_moon import new_system_of_sun_earth_and_moon

def main():
    filename = "../data/SunAndEarthAndMoon_TBBH.amuse"
    ss = new_system_of_sun_earth_and_moon()
    star = ss[0]
    planet = ss[1]
    moon = ss[2]
    converter = nbody_system.nbody_to_si(star.mass, planet.position.length())
    star_gravity = Ph4(converter)
    star_gravity.particles.add_particle(star)

    planet_gravity = Ph4(converter)
    planet_gravity.particles.add_particle(planet)

    moon_gravity = Ph4(converter)
    moon_gravity.particles.add_particle(moon)

    channel_from_star_to_framework = star_gravity.particles.new_channel_to(ss)
    channel_from_planet_to_framework = planet_gravity.particles.new_channel_to(
        ss)
    channel_from_moon_to_framework = moon_gravity.particles.new_channel_to(ss)

    time = 0  |units.yr
    write_set_to_file(ss, filename,
                      timestamp=time, overwrite_file=True)

# #BOOKLISTSTART# #
    sp_gravity = bridge.Bridge()
    sp_gravity.add_system(moon_gravity, (planet_gravity,))
    sp_gravity.add_system(planet_gravity, (moon_gravity,))

    gravity = bridge.Bridge()
    gravity.add_system(sp_gravity, (star_gravity,))
    gravity.add_system(star_gravity, (sp_gravity,))
# #BOOKLISTSTOP# #

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    sp_gravity.timestep = 3 | units.day
    gravity.timestep = 3 | units.day
    dt = 0.1 | units.yr
    t_end = 20 | units.yr
    ddE_max = 0
    while time < t_end:
        time += dt
        gravity.evolve_model(time)

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        channel_from_star_to_framework.copy()
        channel_from_planet_to_framework.copy()
        channel_from_moon_to_framework.copy()
        write_set_to_file(ss, filename, 
                          timestamp=time, append_to_file=True, overwrite_file=True)

        Ekin = gravity.kinetic_energy
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        ddE_max = max(ddE_max, np.abs((Etot_prev-Etot)/Etot))
        print("T=", time, end=' ')
        print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, end=' ')
        print(f"ddE_Max= {ddE_max}")
        Etot_prev = Etot
    gravity.stop()


if __name__ in ('__main__', '__plot__'):
    main()
