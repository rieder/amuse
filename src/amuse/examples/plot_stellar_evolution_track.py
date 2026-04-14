import numpy as np
import matplotlib.pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units, constants
from amuse.community.sse import Sse


def get_stellar_track(mass):

    stars = Particles(1)
    stars[0].mass = mass
    stellar = Sse()
    stellar.particles.add_particles(stars)
    attributes = ["temperature", "luminosity", "stellar_type"]
    to_framework = stellar.particles.new_channel_to(
        stars, attributes=attributes, target_names=attributes
    )
    t = [] | units.Myr
    T = [] | units.K
    L = [] | units.LSun
    stp = [] | units.stellar_type

    Helium_White_Dwarf = 10 | units.stellar_type  # stop here

    while stellar:
        stellar.evolve_model()
        to_framework.copy()
        if stars[0].stellar_type >= Helium_White_Dwarf:
            stellar.stop()
            stellar = False
        else:
            t.append(stellar.model_time)
            T.append(stars[0].temperature)
            L.append(stars[0].luminosity)
            stp.append(stars[0].stellar_type)
            # if mass == 1|units.MSun:
            #     print 'luminosity =', stars[0].luminosity
            #     print 'temperature =', stars[0].temperature
            #     print 'stellar type =', stars[0].stellar_type
            print(stp[-1], stars[0].stellar_type)
            print(stp)

    return t, T, L, stp


sigma = constants.Stefan_hyphen_Boltzmann_constant


def stellar_radius(L, T):
    return np.sqrt(L / (4 * np.pi * sigma * T**4))


def stellar_luminosity(R, T):
    return 4 * np.pi * R**2 * sigma * T**4


def stellar_temperature(R, L):
    return (L / (4 * np.pi * sigma * R**2)) ** (1.0 / 4.0)


def get_color_based_on_stellar_type(istp):
    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    print("Nc=", istp, len(color))
    # return color[istp.value_in(units.stellar_type)]

    if istp.value_in(units.stellar_type) <= 1:
        c = color[0]
    elif istp.value_in(units.stellar_type) <= 2:
        c = color[1]
    elif istp.value_in(units.stellar_type) <= 3:
        c = color[2]
    elif istp.value_in(units.stellar_type) <= 4:
        c = color[3]
    elif istp.value_in(units.stellar_type) <= 5:
        c = color[4]
    elif istp.value_in(units.stellar_type) <= 6:
        c = color[5]
    elif istp.value_in(units.stellar_type) <= 9:
        c = color[5]
    elif istp.value_in(units.stellar_type) <= 12:
        c = color[5]
    elif istp.value_in(units.stellar_type) == 13:
        c = color[5]
    elif istp.value_in(units.stellar_type) == 14:
        c = color[5]
    else:
        c = color[0]
    return c


def main():
    figure = plt.figure()
    ax = figure.add_subplot(111)

    ax.set_xlabel("T [$K$]")
    ax.set_ylabel("L [$L_\odot$]")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim(1e6, 1e3)
    ax.set_ylim(1.0e-1, 3e5)

    L = [0.1, 1.0e5] | units.LSun
    T = stellar_temperature(0.01 | units.RSun, L)
    ax.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c="k")
    ax.text(9e5, 2.0e4, "$0.1 R_\odot$", rotation=-55)

    L = [0.1, 1.0e5] | units.LSun
    T = stellar_temperature(1 | units.RSun, L)
    ax.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c="k")
    ax.text(1.7e4, 50, "$1 R_\odot$", rotation=-55)

    L = [0.1, 1.0e5] | units.LSun
    T = stellar_temperature(100 | units.RSun, L)
    ax.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c="k")
    ax.text(2.85e3, 200, "$100 R_\odot$", rotation=-55)

    t, T, L, stp = get_stellar_track(1 | units.MSun)
    for i in range(len(T) - 2):
        c = get_color_based_on_stellar_type(stp[i])
        ax.plot(
            T[i : i + 2].value_in(units.K),
            L[i : i + 2].value_in(units.LSun),
            lw=4,
            c=c,
            label="$1M_\odot$",
        )

    t, T, L, stp = get_stellar_track(5 | units.MSun)
    ax.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c="k")
    for i in range(len(T) - 2):
        c = get_color_based_on_stellar_type(stp[i])
        ax.plot(
            T[i : i + 2].value_in(units.K),
            L[i : i + 2].value_in(units.LSun),
            lw=4,
            c=c,
            label="$5M_\odot$",
        )

    t, T, L, stp = get_stellar_track(20 | units.MSun)
    ax.plot(T.value_in(units.K), L.value_in(units.LSun), lw=1, c="k")
    for i in range(len(T) - 2):
        c = get_color_based_on_stellar_type(stp[i])
        ax.plot(
            T[i : i + 2].value_in(units.K),
            L[i : i + 2].value_in(units.LSun),
            lw=4,
            c=c,
            label="$20M_\odot$",
        )

    Ttext = 700000
    Ltext = 0.1
    stellar_type = ["", "main sequence", "sub giant", "giant", "remnant"]
    stp = [1, 2, 4, 5, 10] | units.stellar_type
    dL = 2.0
    i = 0
    for sti in range(len(stellar_type)):
        c = get_color_based_on_stellar_type(stp[sti])
        ax.text(Ttext, dL**sti * Ltext, stellar_type[sti], color=c)  # , fontsize=8)

    save_file = "../figures/fig_stellar_evolution_track.pdf"
    plt.savefig(save_file)
    print("\nSaved figure in file", save_file, "\n")
    plt.show()


if __name__ == "__main__":
    main()
