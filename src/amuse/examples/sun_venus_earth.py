"""
Integrate orbits of Venus and the Earth.
"""

import sys
import os
import argparse
import matplotlib.pyplot as plt
from amuse.io import write_set_to_file

# #BOOKLISTSTART1# #
from amuse.datamodel import Particles
from amuse.units import units


def new_sun_venus_and_earth():
    particles = Particles(3)
    sun = particles[0]
    sun.name = "Sun"
    sun.mass = 1.0 | units.MSun
    sun.radius = 1.0 | units.RSun
    sun.position = (855251, -804836, -3186) | units.km
    sun.velocity = (7.893, 11.894, 0.20642) | (units.m / units.s)

    venus = particles[1]
    venus.name = "Venus"
    venus.mass = 0.0025642 | units.MJupiter
    venus.radius = 3026.0 | units.km
    venus.position = (-0.3767, 0.60159, 0.03930) | units.au
    venus.velocity = (-29.7725, -18.849, 0.795) | units.kms

    earth = particles[2]
    earth.name = "Earth"
    earth.mass = 1.0 | units.MEarth
    earth.radius = 1.0 | units.REarth
    earth.position = (-0.98561, 0.0762, -7.847e-5) | units.au
    earth.velocity = (-2.927, -29.803, -0.0005327) | units.kms

    sun.color = "#FFFF00"
    venus.color = "wheat"
    earth.color = "deepskyblue"

    particles.move_to_center()
    return particles


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
from amuse.units import nbody_system
from amuse.community.huayno import Huayno


def integrate_solar_system(
    particles, time_step=1 | units.day, time_end=3 | units.yr,
):
    convert_nbody = nbody_system.nbody_to_si(
        particles.mass.sum(),
        particles[1].position.length(),
    )

    gravity = Huayno(convert_nbody)
    particles_in_code = gravity.particles.add_particles(particles)

    while gravity.model_time < time_end:
        gravity.evolve_model(gravity.model_time + (time_step))
        particles_in_code.new_channel_to(particles).copy()
        particles.savepoint(timestamp=gravity.model_time)

    gravity.stop()

    return particles


# #BOOKLISTSTOP2# #


# #BOOKLISTSTART3# #
def plot_track(particles, output_filename):
    figure = plt.figure(figsize=(10, 10))
    plt.rcParams.update({"font.size": 30})
    ax = figure.add_subplot(1, 1, 1)
    ax.minorticks_on()
    ax.locator_params(nbins=3)

    x_label = "x [au]"
    y_label = "y [au]"
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    ax.scatter([0.0], [0.0], color="#CCCC00", lw=8)
    for planet in particles:
        if planet.name.lower() in ["venus", "earth"]:
            ax.plot(
                planet.get_timeline_of_attribute_as_vector("x")[1].value_in(units.au),
                planet.get_timeline_of_attribute_as_vector("y")[1].value_in(units.au),
                color=planet.color,
            )
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)

    plt.savefig(f"{output_filename}.pdf")
    print(f"\nSaved figure in file {output_filename}\n")


# #BOOKLISTSTOP3# #


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-o",
        "--output_filename",
        default="sun_venus_earth",
        help="output filename",
    )
    result.add_argument(
        "-s",
        "--time_step",
        type=units.day,
        default=1 | units.day,
        help="time step",
    )
    result.add_argument(
        "-t",
        "--time_end",
        type=units.yr,
        default=3.0 | units.yr,
        help="end time of the simulation",
    )
    result.add_argument(
        "-x",
        "--enable_x3d",
        default=False,
        action="store_true",
        help="Create x3d html file",
    )
    return result


def sun_venus_earth(output_filename="sun_venus_earth", enable_x3d=False, **kwargs):
    output_filename_new = output_filename
    i = 0
    while os.path.exists(f"{output_filename_new}.amuse"):
        output_filename_new = f"{output_filename}_{i}"
        i += 1

    particles = new_sun_venus_and_earth()
    particles = integrate_solar_system(particles, **kwargs)

    for savepoint in particles.history:
        write_set_to_file(
            savepoint, f"{output_filename_new}.amuse", append_to_file=True
        )
    plot_track(particles, output_filename)
    if enable_x3d:
        from .plot_x3d import plot_track_x3dom

        plot_track_x3dom(particles, f"{output_filename}.html")


def main(**kwargs):
    sun_venus_earth(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
