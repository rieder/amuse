"""
Minimalistic routine for running the Brutus gravity code
"""

import argparse

import numpy as np
from scipy.optimize import leastsq
from matplotlib import pyplot as plt

from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.brutus.interface import BrutusInterface, Brutus

def phase_space_distance(N, dx, dy, dz, dvx, dvy, dvz):
    return (
        dx.value_in(nbody_system.length) ** 2
        + dy.value_in(nbody_system.length) ** 2
        + dz.value_in(nbody_system.length) ** 2
        + dvx.value_in(nbody_system.speed) ** 2
        + dvy.value_in(nbody_system.speed) ** 2
        + dvz.value_in(nbody_system.speed) ** 2
    ) / (6 * N)


def run_brutus(
    *,
    number_of_particles,
    time_end,
    bs_tolerance,
    word_length,
    dt_param,
    seed,
):

    np.random.seed(seed)
    bodies = new_plummer_model(number_of_particles)
    bodies[0].mass = 1 | nbody_system.mass
    bodies[1].mass = 1 | nbody_system.mass
    bodies[2].mass = 0.1 | nbody_system.mass
    bodies.scale_to_standard()

    # #BOOKLISTSTART# #
    gravity = Brutus()
    gravity.parameters.bs_tolerance = bs_tolerance
    gravity.parameters.word_length = word_length
    gravity.parameters.dt_param = dt_param

    gravity.particles.add_particles(bodies)
    channel_to_framework = gravity.particles.new_channel_to(bodies)
    # #BOOKLISTSTOP# #

    energy_tot_init = gravity.kinetic_energy + gravity.potential_energy

    x = [] | nbody_system.length
    y = [] | nbody_system.length
    z = [] | nbody_system.length
    vx = [] | nbody_system.speed
    vy = [] | nbody_system.speed
    vz = [] | nbody_system.speed
    model_time = 0 | nbody_system.time
    dt = time_end / 1000
    time = [] | nbody_system.time
    while gravity.model_time <= time_end:
        model_time += dt
        gravity.evolve_model(model_time)
        channel_to_framework.copy()
        time.append(model_time)
        x.append(gravity.particles.x)
        y.append(gravity.particles.y)
        z.append(gravity.particles.z)
        vx.append(gravity.particles.vx)
        vy.append(gravity.particles.vy)
        vz.append(gravity.particles.vz)

        energy_kin = gravity.kinetic_energy
        energy_pot = gravity.potential_energy
        energy_tot = energy_kin + energy_pot
        Q = energy_kin / energy_pot
        print(
            f"T= {gravity.model_time} M= {bodies.total_mass()} "
            f"E= {energy_tot} Q= {Q} "
            f"dE= {(energy_tot_init-energy_tot)/energy_tot}"
        )

    gravity.stop()
    return time, x, y, z, vx, vy, vz


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-N",
        default=3,
        dest="number_of_particles",
        type=int,
        help="number of stars",
    )
    result.add_argument(
        "--seed",
        default=3141,
        type=int,
        help="random number seed",
    )
    result.add_argument(
        "-t",
        default=10 | nbody_system.time,
        dest="time_end",
        type=nbody_system.time,
        help="end time of the simulation",
    )
    result.add_argument(
        "--bs_tolerance",
        default=1e-16,
        dest="bs_tolerance",
        type=float,
        help="Bulirsh-Stoer tolerance",
    )
    result.add_argument(
        "-L",
        "--word_length",
        default=180,
        dest="word_length",
        type=float,
        help="Word length",
    )
    result.add_argument(
        "--dt", default=1.0e-5, dest="dt_param", type=float, help="dt parameter"
    )

    return result


def brutus_example(
    number_of_particles=3,
    time_end=1 | nbody_system.time,
    bs_tolerance=1e-30,
    word_length=180,
    dt_param=0.0000000000010,
    seed=3141,
):
    Lw = int(4 * np.abs(np.log10(bs_tolerance)) + 32)
    print(f"Wordlength={Lw}")

    t, x, y, z, vx, vy, vz = run_brutus(
        number_of_particles=number_of_particles,
        time_end=time_end,
        bs_tolerance=bs_tolerance,
        word_length=Lw,
        dt_param=dt_param,
        seed=seed,
    )

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(1, 1, 1)
    ax.scatter(
        x[0].value_in(nbody_system.length),
        y[0].value_in(nbody_system.length),
        s=100,
        c="k",
    )
    ax.plot(x.value_in(nbody_system.length), y.value_in(nbody_system.length), lw=4)
    bs_tolerance = 1.0e-3
    Lw = int(4 * np.abs(np.log10(bs_tolerance)) + 32)
    print(f"Wordlength={Lw}")
    t1, x1, y1, z1, vx1, vy1, vz1 = run_brutus(
        number_of_particles=3,
        time_end=10 | nbody_system.time,
        bs_tolerance=bs_tolerance,
        word_length=Lw,
        dt_param=0.1,
        seed=3141,
    )

    plt.plot(
        x1.value_in(nbody_system.length), y1.value_in(nbody_system.length), lw=1, c="k"
    )
    plt.xlabel("x [N-body units]")
    plt.ylabel("y [N-body units]")
    plt.savefig("fig_brutus_3body_orbit.pdf")
    plt.show()

    fig = plt.figure(figsize=(6, 5))
    dx = x - x1
    dy = y - y1
    dz = z - z1
    plt.plot(dx.value_in(nbody_system.length), dy.value_in(nbody_system.length), lw=4)
    plt.xlabel("dx [N-body units]")
    plt.ylabel("dy [N-body units]")
    plt.show()

    fig = plt.figure(figsize=(6, 5))
    dvx = vx - vx1
    dvy = vy - vy1
    dvz = vz - vz1
    d = phase_space_distance(3, dx, dy, dz, dvx, dvy, dvz)
    time = t.value_in(nbody_system.time)

    plt.plot(time, d[:, 0], lw=4)

    # now fit the power-low to the phase-space distance vs time

    delta = np.log10(d.T[0])
    funcLine = lambda tpl, time: tpl[0] * time + tpl[1]
    func = funcLine
    ErrorFunc = lambda tpl, time, delta: func(tpl, time) - delta
    tplInitial = (1.0, 2.0)
    tplFinal, success = leastsq(ErrorFunc, tplInitial[:], args=(time, delta))
    print("quadratic fit", tplFinal)

    time_fit = np.linspace(time.min(), time.max(), 3)
    psd_fit = 10 ** func(tplFinal, time_fit)
    plt.plot(time_fit, psd_fit, c="k", lw=2, ls="--")
    plt.semilogy()
    plt.xlabel("t [N-body system]")
    plt.ylabel(r"$\log_{10}(\delta)$")
    plt.savefig("fig_brutus_3body_phasespace_distance.pdf")
    plt.show()


def main():
    arguments = new_argument_parser().parse_args()
    brutus_example(**arguments.__dict__)


if __name__ == "__main__":
    main()
