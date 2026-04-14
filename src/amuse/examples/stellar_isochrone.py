"""
Evolve a population of N stars.
initial mass function between Mmin and Mmax and with stellar evolution with
metalicity z.
"""

import argparse
import numpy
import matplotlib.pyplot as plt

from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.seba import Seba
from amuse.community.sse import Sse
from amuse.community.mesa_r2208 import Mesa
from amuse.community.evtwin import Evtwin
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.io import write_set_to_file, read_set_from_file

# from amuse.community.seba.interface import SeBa


def get_stellar_temperature_and_luminosity(
    stars, stellar_code, z=0.02, time_end=100 | units.Myr, write=False, overwrite=False
):
    if stellar_code.lower() in ["seba", "sse"]:
        if stellar_code.lower() == "sse":
            filename = "Stellar_SSE.amuse"
            stellar = Sse()
        else:
            filename = "Stellar_SeBa.amuse"
            stellar = Seba()
        stellar.parameters.metallicity = z
        stellar.particles.add_particles(stars)
        channel_to_framework = stellar.particles.new_channel_to(stars)
        stellar.evolve_model(time_end)
        channel_to_framework.copy_attributes(["radius", "temperature", "luminosity"])
        stellar.stop()
    else:
        if stellar_code.lower() == "mesa":
            filename = "Stellar_MESA.amuse"
        if stellar_code.lower() == "evtwin":
            filename = "Stellar_EVtwin.amuse"

        for si in stars:
            stellar = Mesa()
            stellar.parameters.metallicity = z
            stellar.particles.add_particle(si)
            # stellar.commit_particles()
            channel_to_framework = stellar.particles.new_channel_to(stars)
            try:
                stellar.evolve_model(time_end)
                channel_to_framework.copy_attributes(
                    ["radius", "temperature", "luminosity"]
                )
                print("Successvolly evolved star: m=", si.mass.in_(units.MSun))
            except:
                print("Failed to evolve star: m=", si.mass.in_(units.MSun))
                # stellar.evolve_model(time_end)
            stellar.stop()
    if write:
        write_set_to_file(stars, filename, overwrite_file=overwrite)


def plot_hrd(filename, ax):
    stars = read_set_from_file(filename)
    T = stars.temperature.value_in(units.K)
    L = stars.luminosity.value_in(units.LSun)
    R = stars.radius.value_in(units.RSun)

    R = 80 * numpy.sqrt(R)
    ax.scatter(T, L, lw=0, s=R)


def stellar_isochrone(N, time_end, z, stellar_code="SeBa", plot=False, overwrite=False):
    if plot and "SeBa" not in stellar_code:
        x_label = "T [K]"
        y_label = r"L [L$_\odot$]"
        figure = plt.figure()
        ax = figure.add_subplot(1, 1, 1)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1.0e5, 1.0e3)
        ax.set_ylim(1.0e-4, 1.0e4)
        filename = "Stellar_" + "SeBa" + ".amuse"
        plot_hrd(filename, ax)
        filename = f"Stellar_{stellar_code}.amuse"
        plot_hrd(filename, ax)
        plt.savefig("HRD_N3000at4500Myr")
    elif not plot:
        numpy.random.seed(1)
        masses = new_salpeter_mass_distribution(N)
        stars = Particles(mass=masses)
        get_stellar_temperature_and_luminosity(
            stars,
            stellar_code=stellar_code,
            z=z,
            time_end=time_end,
            write=True,
            overwrite=overwrite,
        )
    else:
        x_label = "T [K]"
        y_label = r"L [L$_\odot$]"
        figure = plt.figure()
        ax = figure.add_subplot(1, 1, 1)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1.0e5, 1.0e3)
        ax.set_ylim(1.0e-4, 1.0e4)
        filename = f"Stellar_{stellar_code}.amuse"
        plot_hrd(filename, ax)
        plt.savefig("HRD_N3000at4500Myr.pdf")


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-C", "--stellar_code", default="SeBa", help="stellar evolution code to use"
    )
    parser.add_argument("-N", type=int, default=3000, help="number of stars")
    parser.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=4500.0 | units.Myr,
        help="end time of the simulation",
    )
    parser.add_argument("-z", type=float, default=0.02, help="metalicity")
    parser.add_argument("-p", "--plot", action="store_true", default=False, help="plot")
    parser.add_argument(
        "--overwrite", action="store_true", default=False, help="overwrite files"
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    stellar_isochrone(**args.__dict__)


if __name__ == "__main__":
    main()
