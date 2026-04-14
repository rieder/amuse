"""
   Script to initialize a star and print its structure
"""

import argparse
import numpy as np
from matplotlib import pyplot as plt

from amuse.units import units
from amuse.datamodel import Particle, Particles
from amuse.community.mesa import Mesa
from amuse.community.mmams import Mmams
from amuse.community.evtwin import Evtwin


def merge_two_stars(Mprim, Msec, t_coll):
    print("Merge two stars.")
    bodies = Particles(
        mass=[Mprim.value_in(units.MSun), Msec.value_in(units.MSun)] | units.MSun
    )

    stellar = Mesa()
    # stellar = EVtwin()
    primary = stellar.particles.add_particles(bodies[0].as_set())
    secondary = stellar.particles.add_particles(bodies[1].as_set())

    print("Evolve for 1Myr")
    stellar.evolve_model(t_coll)

    print("Pre merger:\n", stellar.particles)
    stellar.merge_colliding(
        primary.copy(),
        secondary.copy(),
        Mmams,
        {},
        {"target_n_shells_mixing": 2000},
        return_merge_products=["se"],
    )
    print("Post merger:\n", stellar.particles)

    radius = stellar.particles[0].get_radius_profile()
    rho = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho


def get_density_profile(code=Mesa, M=1.0 | units.MSun, z=0.02):
    stellar = code()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=M))
    print("Nzones=", stellar.particles.get_number_of_zones())
    radius = stellar.particles[0].get_radius_profile()
    rho = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho


def main(M, z, output_filename):
    np.random.seed(31415)
    x_label = r"$R$ [R$_\odot$]"
    y_label = r"$\\rho$ [g/cm$^{3}$]"
    figure = plt.figure(figsize=(8, 8))
    ax = figure.add_subplot(111)
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    plt.xlim(0, 2)
    plt.ylim(1.0e-9, 1.0e2)

    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    r, rho = get_density_profile(Evtwin, M, z)
    plt.plot(
        r.value_in(units.RSun),
        rho.value_in(units.g / units.cm**3),
        label="EVtwin",
        c=color[0],
    )
    r, rho = get_density_profile(Mesa, M, z)
    plt.plot(
        r.value_in(units.RSun),
        rho.value_in(units.g / units.cm**3),
        label="MESA",
        c=color[1],
    )

    # Run the merger code.

    r, rho = merge_two_stars(0.5 * M, 0.5 * M, 1 | units.yr)

    plt.plot(
        r.value_in(units.RSun),
        rho.value_in(units.g / units.cm**3),
        label="MMAMS",
        c=color[2],
    )
    plt.legend(loc="lower right")
    plt.semilogy()

    if output_filename is not None:
        plt.savefig(output_filename)
        print("\nSaved figure in file", output_filename, "\n")
    else:
        output_filename = "../figures/merger_stellar_density_profile.pdf"
        plt.savefig(output_filename)
        print("\nSaved figure in file", output_filename, "\n")
        plt.show()


def new_argument_parser():
    result = argparse.ArgumentParser()
    result.add_option(
        "-M",
        type=units.MSun,
        dest="M",
        default=2.0 | units.MSun,
        help="stellar mass",
    )
    result.add_option(
        "-o", dest="output_filename", default=None, help="output filename"
    )
    result.add_option(
        "-z", dest="z", type=float, default=0.02, help="metallicity"
    )
    return result


if __name__ in ("__main__", "__plot__"):
    args = new_argument_parser().parse_args()
    main(**args.__dict__)
