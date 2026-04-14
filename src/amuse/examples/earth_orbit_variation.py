import argparse
import matplotlib.pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.ext.solarsystem import solar_system_in_time
from amuse.community.huayno import Huayno
from amuse.community.kepler import Kepler


def integrate_solar_system(particles, end_time):
    convert_nbody = nbody_system.nbody_to_si(
        particles.mass.sum(), particles[1].position.length()
    )

    gravity = Huayno(convert_nbody)
    gravity.particles.add_particles(particles)

    SunEarth = Particles()
    SunEarth.add_particle(particles[0])
    SunEarth.add_particle(particles[3])
    channel_from_to_SunEarth = gravity.particles.new_channel_to(SunEarth)

    x_label = "$a-a_0$ [au]"
    y_label = "eccentricty"
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    kep = Kepler(convert_nbody)
    kep.initialize_from_particles(SunEarth)
    a, e = kep.get_elements()
    a0 = a
    ax.scatter((a - a0).value_in(units.au), e, c=colors[1], lw=0, s=20, zorder=20)
    ax.text(
        (a - a0).value_in(units.au) - 0.0003,
        e,
        f"{gravity.model_time.value_in(units.kyr):.0f} kyr",
    )

    dt = 100 | units.yr
    t_diag = 10000 | units.yr
    while gravity.model_time < end_time:

        kep.initialize_from_particles(SunEarth)
        a, e = kep.get_elements()
        ax.scatter((a - a0).value_in(units.au), e, c=colors[0], lw=0, s=20, zorder=1)

        time = gravity.model_time
        if time > t_diag:
            ax.scatter((a - a0).value_in(units.au), e, c=colors[1], lw=0, s=20, zorder=20)
            t_diag += 10000 | units.yr
            fmt = f"{time.value_in(units.kyr):.0f} kyr"
            if time < 100 | units.yr:
                ax.text(-0.0006, e, fmt)
            elif time < 15000 | units.yr:
                ax.text(0.0003, e, fmt)
            elif time < 25000 | units.yr:
                ax.text(0.0017, e, fmt)
            elif time < 35000 | units.yr:
                ax.text(0.0031, e, fmt)
            elif time < 45000 | units.yr:
                ax.text(0.0028, e, fmt)
            elif time < 55000 | units.yr:
                ax.text(0.0021, e, fmt)
            elif time < 65000 | units.yr:
                ax.text(0.0017, e + 0.002, fmt)
            elif time < 75000 | units.yr:
                ax.text(0.0014, e - 0.002, fmt)
            else:
                pass

        gravity.evolve_model(gravity.model_time + dt)
        channel_from_to_SunEarth.copy()
    gravity.stop()

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim((-0.001, 0.004))

    file = "EarthOrbitVariation.png"
    plt.savefig(file)
    print("\nSaved figure in file", file, "\n")

    return


def earth_orbit_variation(epoch=2474649.5 | units.day, time_end=80000 | units.yr):
    particles = solar_system_in_time(time_JD=epoch)
    integrate_solar_system(particles, time_end)
    return


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-t",
        "--time_end",
        type=units.yr,
        default=80000 | units.yr,
        help="end time",
    )
    parser.add_argument(
        "-e",
        "--epoch",
        type=units.day,
        default=2474649.5 | units.day,
        help="epoch",
    )
    return parser


def main(**kwargs):
    earth_orbit_variation(**kwargs)


if __name__ == "__main__":
    main()
