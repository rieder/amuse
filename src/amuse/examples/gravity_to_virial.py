"""
Simple routine for running a gravity code
"""

import argparse
import matplotlib.pyplot as plt
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.community.huayno import Huayno
from amuse.community.ph4 import Ph4
from amuse.community.bhtree import Bhtree


def virial_ratio_evolution(code, bodies, Q_init, time_end):
    dt = 0.06125 | time_end.unit
    bodies.scale_to_standard(virial_ratio=Q_init)
    bodies.radius = 0 | nbody_system.length
    gravity = code()
    gravity.particles.add_particles(bodies)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)

    Etot_prev = Etot_init = gravity.kinetic_energy + gravity.potential_energy
    time = [0.0] | time_end.unit
    Q = [Q_init]
    while time[-1] < time_end:
        time.append(time[-1] + dt)
        gravity.evolve_model(time[-1])
        channel_from_gravity_to_framework.copy()
        Ekin = gravity.kinetic_energy
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        Q.append(-1 * Ekin / Epot)
        print("T=", time[-1], "Q= ", Q[-1], end=" ")
        print("M=", bodies.mass.sum(), "E= ", Etot, end=" ")
        print("dE=", (Etot_init - Etot) / Etot, "ddE=", (Etot_prev - Etot) / Etot)
        Etot_prev = Etot
    gravity.stop()
    return time, Q


def gravity_to_virial(number_of_stars, time_end):
    Q_init = 0.2
    particles = new_plummer_model(number_of_stars)
    codes = [Ph4, Huayno, Bhtree]
    ci = 0
    x_label = "time [N-body units]"
    # y_label = "$Q [\equiv -E_{\rm kin}/E_{\rm pot}]$"
    y_label = "virial ratio $Q$"
    figure = plt.figure()
    ax1 = figure.add_subplot(1, 1, 1)
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)
    ax1.set_xlim(0, time_end.value_in(time_end.unit))
    ax1.set_ylim(0, 0.65)
    ax1.plot([0, time_end.value_in(time_end.unit)], [0.5, 0.5], lw=1, ls="--", c="k")
    for code in codes:
        time, Q = virial_ratio_evolution(code, particles, Q_init, time_end)
        ax1.plot(
            time.value_in(time_end.unit),
            Q,
        )
        ci += 1
    plt.savefig("gravity_to_virial")


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-N",
        "--number_of_stars",
        type=int,
        default=1000,
        help="number of stars",
    )
    parser.add_argument(
        "-t",
        "--time_end",
        type=nbody_system.time,
        default=2 | nbody_system.time,
        help="end time of the simulation",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    gravity_to_virial(**args.__dict__)


if __name__ == "__main__":
    main()
