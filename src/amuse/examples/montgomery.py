import numpy as np
import matplotlib.pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units, constants, nbody_system
from amuse.community.hermite_grx import Hermitegrx


def get_initial_conditions_figure_eight(scale_mass, scale_length):
    unit_mass = scale_mass.as_unit()
    unit_length = scale_length.as_unit()
    converter = nbody_system.nbody_to_si(scale_mass, scale_length)

    particles = Particles(3)
    particles.position = (
        np.array(
            [
                [0.9700436, -0.24308753, 0],
                [-0.9700436, 0.24308753, 0],
                [0, 0, 0],
            ]
        )
        | unit_length
    )
    particles.velocity = np.array(
        [
            [0.466203685, 0.43236573, 0],
            [0.466203685, 0.43236573, 0],
            [-2 * 0.466203685, -2 * 0.43236573, 0],
        ]
    ) | (unit_length / converter.to_si(nbody_system.time).as_unit())
    particles.mass = np.array([1, 1, 1]) | unit_mass
    particles.radius = 0 | unit_length

    particles.move_to_center()

    return particles, converter


def get_energy(grav, pert=None):
    if pert is None:
        return grav.particles.kinetic_energy() + grav.particles.potential_energy()
    return grav.get_total_energy_with(pert)[0]


def get_trajectories(initial, grav, time_end, pert=None, number_of_steps=10000):
    grav.particles.add_particles(initial[0])

    # pre_energy = get_energy(grav, pert)
    x1 = [] | units.au
    x2 = [] | units.au
    x3 = [] | units.au
    y1 = [] | units.au
    y2 = [] | units.au
    y3 = [] | units.au
    z1 = [] | units.au
    z2 = [] | units.au
    z3 = [] | units.au
    times = [] | units.yr
    energies = [] | units.J

    time = 0 * time_end
    time_step = time_end / number_of_steps

    while time < time_end:
        x1.append(grav.particles[0].position.x)
        y1.append(grav.particles[0].position.y)
        z1.append(grav.particles[0].position.z)
        x2.append(grav.particles[1].position.x)
        y2.append(grav.particles[1].position.y)
        z2.append(grav.particles[1].position.z)
        x3.append(grav.particles[2].position.x)
        y3.append(grav.particles[2].position.y)
        z3.append(grav.particles[2].position.z)

        times.append(time)
        energies.append(get_energy(grav, pert))

        time += time_step
        # print("Integrate:", time)
        grav.evolve_model(time)

    # grav.close()

    return x1, y1, z1, x2, y2, z2, x3, y3, z3, times, energies


def run_montgomery(relative_light_speed=1):
    initial = get_initial_conditions_figure_eight(
        1 | units.MSun,
        1 | units.RSun,
    )
    # bodies = initial[0]
    converter = initial[1]
    # #BOOKLISTSTART1# #
    n_body_code = Hermitegrx
    grav = n_body_code(converter)
    pert = "2.5PN_EIH"
    grav.parameters.perturbation = pert
    grav.parameters.integrator = "RegularizedHermite"
    grav.parameters.dt_param = 0.01
    grav.parameters.light_speed = relative_light_speed * constants.c
    # #BOOKLISTSTOP1# #

    x1, y1, z1, x2, y2, z2, x3, y3, z3, times, energies = get_trajectories(
        initial, grav, 0.001 | units.yr, pert
    )
    grav.stop()
    return x1, y1, z1, x2, y2, z2, x3, y3, z3


def montgomery():
    unit_length = units.au
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)

    x1, y1, z1, x2, y2, z2, x3, y3, z3 = run_montgomery(1.0)
    ax.plot(x1.value_in(unit_length), y1.value_in(unit_length), c=color[0], lw=6)
    ax.plot(x2.value_in(unit_length), y2.value_in(unit_length), c=color[1], lw=6)
    ax.plot(x3.value_in(unit_length), y3.value_in(unit_length), c=color[2], lw=6)

    x1, y1, z1, x2, y2, z2, x3, y3, z3 = run_montgomery(0.01)
    ax.plot(x1.value_in(unit_length), y1.value_in(unit_length), c=color[0], lw=1)
    ax.plot(x2.value_in(unit_length), y2.value_in(unit_length), c=color[1], lw=1)
    ax.plot(x3.value_in(unit_length), y3.value_in(unit_length), c=color[2], lw=1)

    ax.set_xlabel(r"$x [R_\odot]$")
    ax.set_ylabel(r"$y [R_\odot]$")

    plt.savefig("fig_Montgomery_GRX.pdf")


def main(**kwargs):
    montgomery(**kwargs)


if __name__ == "__main__":
    main()
