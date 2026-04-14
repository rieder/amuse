import argparse
from matplotlib import pyplot as plt
from amuse.community.mesa import Mesa
from amuse.units import units
from amuse.datamodel import Particle

import pickle

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type
HeWhiteDwarf = 10 | units.stellar_type


def stellar_remnant(stellar):
    remnant = True
    if (
        stellar.particles[0].stellar_type < HeWhiteDwarf
        or stellar.particles[0].stellar_type > 11 | units.stellar_type
    ):
        remnant = False
    return remnant


def stellar_core_temperature_and_density(mass, metallicity=0.02):
    filename = f"stellar_M{mass.value_in(units.MSun)}MSun_core_temperature_and_density.pkl"

    time = [] | units.Myr
    rho_core = [] | units.g / units.cm**3
    T_core = [] | units.K
    stellar_type = [] | units.stellar_type
    check = 0
    ncheck = 100

    stellar = Mesa(version="15140")
    stellar.parameters.metallicity = metallicity
    star = stellar.particles.add_particle(Particle(mass=mass))

    # #BOOKLISTSTART1# #
    while not stellar_remnant(stellar):

        star.evolve_one_step()

        nzones = star.get_number_of_zones()
        rhoc = star.get_density_profile(nzones)[0]
        Tc = star.get_temperature_profile(nzones)[0]

        time.append(stellar.model_time)
        rho_core.append(rhoc)
        T_core.append(Tc)
        stellar_type.append(star.stellar_type)
        si = star.stellar_type.value_in(units.stellar_type)

        with open(filename, "wb") as file:
            pickle.dump([time, rho_core, T_core, stellar_type], file)
        # #BOOKLISTSTOP1# #

        check += 1
        if check == ncheck:
            check = 0
            try:
                x = open("STOP")
                stop = True
            except:
                stop = False
            if stop:
                break

        # if time[-1] >= 3.1339239457|units.Myr: break
        print(
            star.age.in_(units.Myr),
            rhoc.in_(units.g / units.cm**3),
            Tc.in_(units.K),
            star.stellar_type,
        )

    stellar.stop()

    return time, rho_core, T_core, stellar_type


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-m",
        "--mass",
        type=units.MSun,
        default=1 | units.MSun,
        help="Stellar mass",
    )
    parser.add_argument(
        "-Z",
        "--metallicity",
        default=0.02,
        help="Stellar metallicity",
    )
    return parser


def calculate_core_temperature_density(mass=1 | units.MSun, metallicity=0.02):
    filename = (
        f"stellar_M{mass.value_in(units.MSun)}MSun_core_temperature_and_density.pkl"
    )

    time, rhoc, Tc, stp = stellar_core_temperature_and_density(mass, metallicity)
    with open(filename, "wb") as file:
        pickle.dump([time, rhoc, Tc, stp], file)

    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.scatter(rhoc.value_in(units.g / units.cm**3), Tc.value_in(units.K))

    fontsize = 12
    ax.text(20.0, 1.5e7, r"$1\,M_\odot$", fontsize=fontsize)
    ax.text(1.0e2, 7.0e7, r"$10\,M_\odot$", fontsize=fontsize)
    ax.text(4.0e4, 1.6e9, r"$100\,M_\odot$", fontsize=fontsize)

    ax.set_xlabel("core density [g/cm$^3$]")
    ax.set_ylabel("core temperature [K]")
    ax.set_xlim((1.0, 1.0e7))
    ax.set_ylim((1.0e7, 1.0e10))
    ax.set_xscale("log")
    ax.set_yscale("log")

    save_file = "calculate_core_temperature_density.pdf"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}")


def main(**kwargs):
    calculate_core_temperature_density(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
