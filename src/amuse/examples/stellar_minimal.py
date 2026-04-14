"""
Minimal example for running a stellar evolution code in AMUSE
"""

import sys
import argparse
import numpy as np

# #BOOKLISTSTART# #
from amuse.units import units
from amuse.datamodel import Particle
from amuse.community.seba import Seba


def stellar_minimal(
    mass=1.0 | units.MSun, metallicity=0.02, time_end=4700 | units.Myr
):
    stellar = Seba()
    stellar.parameters.metallicity = metallicity
    star = stellar.particles.add_particle(Particle(mass=mass))

    initial_luminosity = star.luminosity
    time_step = 0.1 | units.Myr

    while stellar.model_time < time_end:
        stellar.evolve_model(stellar.model_time + time_step)

        print(
            f"at T={stellar.model_time.in_(units.Myr)} "
            f"L(t=0)={initial_luminosity}, "
            f"L (t={star.age.in_(units.Myr)})={star.luminosity.in_(units.LSun)}, "
            f"m={star.mass.in_(units.MSun)}, "
            f"R={star.radius.in_(units.RSun)}"
        )

    luminosity = star.luminosity

    print(
        f"Time={stellar.model_time.in_(units.Myr)} "
        f"L={np.log10(luminosity.value_in(units.erg / units.s))}"
    )

    stellar.stop()


# #BOOKLISTSTOP# #


# #BOOKLISTSTART2# #
def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-m",
        type=units.MSun,
        dest="mass",
        default=1.0 | units.MSun,
        help="stellar mass",
    )
    parser.add_argument(
        "-t",
        type=units.Myr,
        dest="time_end",
        default=4700.0 | units.Myr,
        help="end time of the simulation",
    )
    parser.add_argument(
        "-z",
        dest="metallicity",
        type=float,
        default=0.02,
        help="metallicity",
    )
    return parser


def main(**kwargs):
    stellar_minimal(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
# #BOOKLISTSTOP2# #
