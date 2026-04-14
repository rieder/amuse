import argparse
import matplotlib.pyplot as plt

from amuse.datamodel import Particle
from amuse.plot import plot
from amuse.community.evtwin import Evtwin
from amuse.community.mesa import Mesa
from amuse.units import units


# #BOOKLISTSTART# #
def get_density_profile(code=Evtwin, mass=2.0 | units.MSun, metallicity=0.02):
    stellar = code()
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particle(Particle(mass=mass))
    radius = stellar.particles[0].get_radius_profile()
    rho = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho


# #BOOKLISTSTOP# #


def initialize_single_star(mass=1.0 | units.MSun, metallicity=0.02):
    r, rho = get_density_profile(Evtwin, mass, metallicity)
    plot(r.in_(units.RSun), rho, label="EVtwin")
    r, rho = get_density_profile(Mesa, mass, metallicity)
    plot(r.in_(units.RSun), rho, label="MESA")
    ax = plt.gca()
    ax.set_xlabel(r"$R$ [$R_\odot$]")
    ax.set_ylabel("density [$g/cm^3$]")
    plt.semilogy()

    save_file = "initialize_single_star.png"
    plt.savefig(save_file)
    print(f"\nSaved figure in file {save_file}\n")


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
        "-z", "--metallicity", type=float, default=0.02, help="metallicity"
    )
    return result


def main(**kwargs):
    initialize_single_star(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
