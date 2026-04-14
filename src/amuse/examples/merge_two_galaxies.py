import argparse
import numpy as np
import matplotlib.pyplot as plt
from amuse.ext.galactics_model import new_galactics_model
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.community.gadget2 import Gadget2


def make_plot(disk1, disk2, filename):
    x_label = "X [kpc]"
    y_label = "Y [kpc]"
    figure = plt.figure(figsize=(8, 8))
    ax = figure.add_subplot(111)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    c = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)

    plt.scatter(
        disk1.x.value_in(units.kpc),
        disk1.y.value_in(units.kpc),
        c=c[0],
        alpha=1,
        s=1,
        lw=0,
    )
    plt.scatter(
        disk2.x.value_in(units.kpc),
        disk2.y.value_in(units.kpc),
        c=c[1],
        alpha=1,
        s=1,
        lw=0,
    )
    plt.savefig(filename)
    print(f"Saved figure in file {filename}")


def make_galaxies(M_galaxy, R_galaxy, n_halo, n_bulge, n_disk):
    converter = nbody_system.nbody_to_si(M_galaxy, R_galaxy)
    galaxy1 = new_galactics_model(
        n_halo,
        converter,
        # do_scale = True,
        bulge_number_of_particles=n_bulge,
        disk_number_of_particles=n_disk,
    )
    galaxy2 = Particles(len(galaxy1))
    galaxy2.mass = galaxy1.mass
    galaxy2.position = galaxy1.position
    galaxy2.velocity = galaxy1.velocity

    galaxy1.rotate(0.0, np.pi / 2, np.pi / 4)
    galaxy1.position += [100.0, 100, 0] | units.kpc
    #    galaxy1.velocity += [-3000.0, 0.0, -3000.0] | units.km/units.s
    galaxy1.velocity += [-10.0, 0.0, -10.0] | units.km / units.s

    galaxy2.rotate(np.pi / 4, np.pi / 4, 0.0)
    galaxy2.position -= [100.0, 0, 0] | units.kpc
    galaxy2.velocity -= [0.0, 0.0, 0] | units.km / units.s

    return galaxy1, galaxy2, converter


def simulate_merger(galaxy1, galaxy2, converter, n_halo, t_end):
    converter = nbody_system.nbody_to_si(1.0e12 | units.MSun, 100 | units.kpc)
    dynamics = Gadget2(converter, number_of_workers=4)
    dynamics.parameters.epsilon_squared = (100 | units.parsec) ** 2
    set1 = dynamics.particles.add_particles(galaxy1)
    set2 = dynamics.particles.add_particles(galaxy2)
    dynamics.particles.move_to_center()
    disk1 = set1[:n_halo]
    disk2 = set2[:n_halo]

    make_plot(disk1, disk2, "Galaxy_merger_t0Myr")
    dynamics.evolve_model(t_end)
    make_plot(disk1, disk2, "Galaxy_merger_t" + str(t_end.value_in(units.Myr)) + "Myr")

    dynamics.stop()


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-M",
        "--mass_galaxy",
        type=units.MSun,
        default=1.0e12 | units.MSun,
        help="Galaxy mass",
    )
    parser.add_argument(
        "-R",
        "--radius_galaxy",
        type=units.kpc,
        default=10 | units.kpc,
        help="Galaxy size",
    )
    parser.add_argument(
        "--n_bulge",
        default=10000,
        help="number of stars in the bulge",
    )
    parser.add_argument(
        "--n_disk",
        default=10000,
        help="number of stars in the disk",
    )
    parser.add_argument(
        "--n_halo",
        default=20000,
        help="number of stars in the halo",
    )
    parser.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=200 | units.Myr,
        help="End of the simulation",
    )
    return parser


def merge_two_galaxies(
    mass_galaxy=1.0e12 | units.MSun,
    radius_galaxy=10 | units.kpc,
    n_bulge=10000,
    n_disk=10000,
    n_halo=20000,
    time_end=200 | units.Myr,
):
    galaxy1, galaxy2, converter = make_galaxies(
        mass_galaxy, radius_galaxy, n_halo, n_bulge, n_disk
    )
    simulate_merger(galaxy1, galaxy2, converter, n_halo, time_end)


def main(**kwargs):
    merge_two_galaxies(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
