import matplotlib.pyplot as plt

from amuse.plot import xlabel, ylabel, effective_iso_potential_plot
from amuse.units import units, constants, nbody_system
from amuse.community.hermite import Hermite
from amuse.datamodel import Particles


def new_sun_earth_system():
    particles = Particles(2)
    particles.mass = [1, 0.2] | units.MSun
    particles.position = [[0, 0, 0], [1.0, 0, 0]] | units.au
    particles.velocity = [0, 0, 0] | units.kms
    particles[1].vy = (constants.G * particles.total_mass() / particles[1].x).sqrt()
    return particles


def setup_gravity_code():
    converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.au)
    gravity = Hermite(converter)
    gravity.particles.add_particles(new_sun_earth_system())
    return gravity


def make_effective_iso_potential_plot(gravity_code):
    omega = (
        constants.G * gravity_code.particles.total_mass() / (1.0 | units.au**3)
    ).sqrt()
    center_of_mass = gravity_code.particles.center_of_mass()[:2]

    figure = plt.figure(figsize=(8, 8))
    ax = figure.add_subplot(1, 1, 1)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    ax.minorticks_on()

    current_axes = plt.subplot(1, 1, 1)
    current_axes.set_aspect("equal", adjustable="box")
    lim = 1.7
    effective_iso_potential_plot(
        gravity_code,
        omega,
        xlim=[-lim, lim] | units.au,
        ylim=[-lim, lim] | units.au,
        center_of_rotation=center_of_mass,
        fraction_screen_filled=0.8,
    )
    xlabel("x")
    ylabel("y")
    plt.text(0.6, -0.06, "$L_1$")
    plt.text(1.35, -0.06, "$L_2$")
    plt.text(-0.99, -0.06, "$L_3$")
    plt.text(0.40, 0.82, "$L_4$")
    plt.text(0.40, -0.92, "$L_5$")

    save_file = "../figures/lagrange_points.pdf"
    plt.savefig(save_file)
    print("\nOutput saved in", save_file)
    plt.show()


def lagrange_points():
    gravity = setup_gravity_code()
    make_effective_iso_potential_plot(gravity)
    gravity.stop()


def main():
    lagrange_points()


if __name__ == "__main__":
    main()
