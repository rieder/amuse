"""
Calculate the response of a star as a result of mass loss.
"""

import sys
import argparse
from amuse.units import units
from amuse.datamodel import Particles
from amuse.support.console import set_printing_strategy
from amuse.community.mesa import Mesa

Second_Asymptotic_Giant_Branch = 6 | units.stellar_type

set_printing_strategy(
    "custom",
    preferred_units=[units.MSun, units.RSun, units.Myr, units.MSun / units.yr],
    precision=6,
    prefix="",
    separator=" [",
    suffix="]",
)


# #BOOKLISTSTART1# #
def calculate_zeta(star, metallicity, dmdt):
    stellar = Mesa()
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particles(star)
    radius_old = star.radius
    star.mass_change = dmdt
    mass_step = 0.01 * star.mass
    star.time_step = mass_step / dmdt
    stellar.particles.evolve_one_step()
    radius_new = stellar.particles[0].radius
    dlnr = (radius_new - radius_old) / radius_old
    dlnm = (stellar.particles[0].mass - star.mass) / star.mass
    zeta = dlnr / dlnm
    stellar.stop()
    return zeta


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
def massloss_response(
    mass=1.0 | units.MSun, metallicity=0.02, dmdt=-0.01 | units.MSun / units.yr
):
    stellar = Mesa()
    stellar.parameters.metallicity = metallicity
    bodies = Particles(mass=mass)
    stellar.particles.add_particles(bodies)
    stellar = turnon_massloss(stellar)
    channel_to_framework = stellar.particles.new_channel_to(bodies)
    copy_argument = ["age", "mass", "radius", "stellar_type"]
    while stellar.particles[0].stellar_type < Second_Asymptotic_Giant_Branch:
        stellar.particles.evolve_one_step()
        channel_to_framework.copy_attributes(copy_argument)
        star = stellar.particles.copy()
        zeta = calculate_zeta(star, metallicity, dmdt)
        print(
            "Zeta=",
            zeta[0],
            bodies[0].age,
            bodies[0].mass,
            bodies[0].radius,
            dmdt,
            bodies[0].stellar_type,
        )
    stellar.stop()


# #BOOKLISTSTOP2# #


def turnon_massloss(stellar):
    if stellar.mesa_version == "2208":
        return stellar  # Mass loss defaults to on in this version
    for particle in stellar.particles:
        particle.set_control("cool_wind_RGB_scheme", "Reimers")
        particle.set_control("Reimers_scaling_factor", 0.1)
        particle.set_control("cool_wind_AGB_scheme", "Blocker")
        particle.set_control("Blocker_scaling_factor", 0.5)
        particle.set_control("RGB_to_AGB_wind_switch", 10**-4)
    return stellar


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-M",
        "--mass",
        type=units.MSun,
        default=1.0 | units.MSun,
        help="stellar mass",
    )
    result.add_argument(
        "--dmdt",
        type=units.MSun / units.yr,
        dest="dmdt",
        default=-0.01 | (units.MSun / units.yr),
        help="dmdt",
    )
    result.add_argument(
        "-z", dest="metallicity", type=float, default=0.02, help="metallicity"
    )
    return result


def main(**kwargs):
    massloss_response(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
