import argparse
import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.seba import Seba
from amuse.community.sse import Sse
from amuse.community.mesa import Mesa
from amuse.community.evtwin import Evtwin


def plot_solar_comparison(time_end, mass, z, Tstar, Lstar):
    stellar_evolution_codes = ["SeBa", "SSE", "MESA", "EVtwin"]
    marker = ["o", "v", "<", ">"]

    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    plt.xlabel("$(T-T_\odot)/T_\odot)$")
    plt.ylabel("$(L-L_\odot)/L_\odot)$")
    plt.xlim(-0.006, 0.004)
    plt.ylim(-0.1, 0.1)
    plt.scatter([0], [0], marker="o", label="Sun", s=200, lw=0)

    for si, code in enumerate(stellar_evolution_codes):
        if code == "SeBa":
            stellar = Seba()
        elif code == "SSE":
            stellar = Sse()
        elif code == "MESA":
            stellar = Mesa()
        elif code == "EVtwin":
            stellar = Evtwin()
        else:
            raise ValueError(f"Unknown stellar evolution code {code}")
        stellar.parameters.metallicity = z

        star = Particles(1)
        star.mass = mass
        # time_end = 6000.0 | units.Myr
        stellar.particles.add_particles(star)
        attributes = ["temperature", "luminosity", "age"]
        to_framework = stellar.particles.new_channel_to(
            star, attributes=attributes, target_names=attributes
        )
        t = [] | units.Myr
        T = []
        L = []
        min_dist_sun = 10000.0
        current_time = 3000.0 | units.Myr

        print(code)
        # dt = 50 | units.Myr
        # time = 4000 | units.Myr
        while stellar:
            print(stellar.model_time.value_in(units.Myr))
            current_time = current_time + stellar.particles[0].time_step
            stellar.evolve_model(current_time)

            to_framework.copy()

            if star[0].age >= time_end:
                stellar.stop()
                stellar = False
            else:
                L.append((star[0].luminosity - Lstar) / Lstar)
                T.append((star[0].temperature - Tstar) / Tstar)
                t.append(star[0].age)

                deltaL = np.abs((star[0].luminosity - Lstar) / Lstar)
                deltaT = np.abs((star[0].temperature - Tstar) / Tstar)

                dist = np.sqrt(deltaL * deltaL + deltaT * deltaT)

                if min_dist_sun > dist:

                    min_dist_sun = dist

                    L_sim_sun = (star[0].luminosity - Lstar) / Lstar
                    T_sim_sun = (star[0].temperature - Tstar) / Tstar

                    eta = star[0].age
        print(eta)
        if si == 3:
            plt.plot(T, L, ls="-", marker=marker[si], markersize=10)
            plt.scatter(
                T_sim_sun, L_sim_sun, marker=marker[si], label=code, s=300, lw=1
            )
        else:
            plt.plot(T, L, ls="-", marker=marker[si], markersize=10)
            plt.scatter(
                T_sim_sun, L_sim_sun, marker=marker[si], label=code, s=300, lw=1
            )

    plt.legend(scatterpoints=1, loc="best")

    save_file = "fig_SunComparison.png"
    plt.savefig(save_file)
    print("\nSaved figure in file", save_file, "\n")
    plt.show()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-T",
        "--Tstar",
        type=units.K,
        default=5778 | units.K,
        help="stellar temperature",
    )
    parser.add_argument(
        "-L",
        "--Lstar",
        type=units.LSun,
        default=1 | units.LSun,
        help="stellar luminosity",
    )
    parser.add_argument(
        "-m", "--mass", type=units.MSun, default=1.0 | units.MSun, help="stellar mass"
    )
    parser.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=4.6 | units.Gyr,
        help="end time of the simulation",
    )
    parser.add_argument("-z", type=float, default=0.02, help="metallicity")
    return parser


def main():
    args = new_argument_parser().parse_args()
    plot_solar_comparison(**args.__dict__)


if __name__ == "__main__":
    main()
