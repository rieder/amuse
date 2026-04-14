import argparse
import numpy as np
import matplotlib.pyplot as plt

from amuse.units import units, nbody_system
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ic.fractalcluster import new_fractal_cluster_model
from amuse.ic.kroupa import new_kroupa_mass_distribution

from amuse.community.hop import Hop


def plot_single_image(groups_of_particles, lim=10):
    left, width = 0.1, 0.4
    bottom, height = 0.1, 0.4
    bottom_h = left_h = left + width + 0.05
    rect_xy = [left, bottom, width, height]
    rect_xz = [left, bottom_h, width, 0.4]
    rect_yz = [left_h, bottom, 0.4, height]

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig = plt.figure(figsize=(10, 10))

    xy = plt.axes(rect_xy)
    xz = plt.axes(rect_xz)
    yz = plt.axes(rect_yz)
    xy.set_xlabel("X [pc]")
    xy.set_ylabel("Y [pc]")
    xz.set_ylabel("Z [pc]")
    yz.set_xlabel("Z [pc]")

    i = 0
    for group in groups_of_particles:
        x = group.x.value_in(units.parsec)
        y = group.y.value_in(units.parsec)
        z = group.z.value_in(units.parsec)
        xy.scatter(x, y, lw=0, c=colors[min(9, i)], s=8)
        xz.scatter(x, z, lw=0, c=colors[min(9, i)], s=8)
        yz.scatter(z, y, lw=0, c=colors[min(9, i)], s=8)
        i += 1

    xy.set_xlim((-lim, lim))
    xy.set_ylim((-lim, lim))
    xz.set_xlim(xy.get_xlim())
    yz.set_ylim(xy.get_xlim())
    yz.set_xlim(xy.get_xlim())
    xz.set_ylim(xy.get_xlim())

    save_file = "FractalClusterHop.pdf"
    plt.savefig(save_file)
    print("\nSaved figure in file", save_file, "\n")
    plt.show()


def find_clumps_with_hop(particles, unit_converter):
    # #BOOKLISTSTART# #
    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()

    mean_density = hop.particles.density.mean()
    hop.parameters.peak_density_threshold = mean_density
    hop.parameters.saddle_density_threshold = 0.99 * mean_density
    hop.parameters.outer_density_threshold = 0.01 * mean_density

    hop.do_hop()
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    hop.stop()
    # #BOOKLISTSTOP# #

    return result


def plot_fractal_clumpy_cluster(
    number_of_particles=2000,
    virial_radius=0.5 | units.parsec,
    virial_ratio=0.5,
    fractal_dimension=1.6,
    seed=12345,
):
    np.random.seed(seed)

    masses = new_kroupa_mass_distribution(number_of_particles)

    converter = nbody_system.nbody_to_si(masses.sum(), virial_radius)
    bodies = new_fractal_cluster_model(
        N=number_of_particles,
        fractal_dimension=fractal_dimension,
        random_seed=seed,
        convert_nbody=converter,
    )
    bodies.mass = masses
    bodies.move_to_center()
    bodies.scale_to_standard(converter, virial_ratio=virial_ratio)

    clumps = find_clumps_with_hop(bodies, converter)

    for clump in clumps:
        clump.scale_to_standard(converter, virial_ratio=0.5)

    plot_single_image(clumps, lim=10)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-N", "--number_of_particles", type=int, default=2000, help="number of stars"
    )
    parser.add_argument(
        "-R",
        "--virial_radius",
        type=units.parsec,
        default=0.5 | units.parsec,
        help="cluster virial radius",
    )
    parser.add_argument(
        "-Q", "--virial_ratio", type=float, default=0.5, help="virial ratio"
    )
    parser.add_argument(
        "-F", "--fractal_dimension", type=float, default=1.6, help="fractal dimension"
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=12345,
        help="random number seed",
    )
    return parser


def main(**kwargs):
    plot_fractal_clumpy_cluster(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
