"""
Simple script to run a `2+2' quadruple system with SecularMultiple.
The two `inner' binaries are denoted with `A' and `B'; the wider outer binary
('superorbit') is denoted with `C'.
Orbital parameters can be provided with command line arguments.
Note: setting N_output to a large value will slow down the script due to Python
overhead, but will make nicer-looking plots.

Adrian Hamers, December 2017
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from amuse.community.secularmultiple import Secularmultiple
from amuse.units import units
from amuse.units.trigo import sin, cos, arccos
from amuse.units.quantities import as_vector_quantity
from amuse.datamodel import Particles


# #BOOKLISTSTART1# #
def initialize_multiple_system(
    number_of_bodies,
    masses,
    semimajor_axis,
    eccentricity,
    inclination,
    argument_of_pericenter,
    longitude_of_ascending_node,
):
    """
    Initializes a system of multiples
    """

    number_of_binaries = number_of_bodies - 1
    particles = Particles(number_of_bodies + number_of_binaries)
    for index in range(number_of_bodies):
        particle = particles[index]
        particle.mass = masses[index]
        particle.is_binary = False
        particle.radius = 1.0 | units.RSun
        particle.child1 = None
        particle.child2 = None

    for index in range(number_of_binaries):
        particle = particles[index + number_of_bodies]
        particle.is_binary = True
        particle.semimajor_axis = semimajor_axis[index]
        particle.eccentricity = eccentricity[index]
        particle.inclination = inclination[index]
        particle.argument_of_pericenter = argument_of_pericenter[index]
        particle.longitude_of_ascending_node = longitude_of_ascending_node[index]

        # Specify the `2+2' hierarchy:

        if index == 0:
            particle.child1 = particles[0]
            particle.child2 = particles[1]
        elif index == 1:
            particle.child1 = particles[2]
            particle.child2 = particles[3]
        elif index == 2:
            particle.child1 = particles[4]
            particle.child2 = particles[5]
    binaries = particles[particles.is_binary]

    return particles, binaries


# #BOOKLISTSTOP1# #


def evolve_quadruple(
    number_of_output_steps=400,
    time_end=5.0 | units.Myr,
    masses=[1.0, 0.8, 1.1, 0.9] | units.MSun,
    semimajor_axis=[1.0, 1.2, 100.0] | units.au,
    eccentricity=[0.1, 0.1, 0.3],
    inclination=[75.0, 80.0, 0.001] | units.deg,
    argument_of_pericenter=[10.0, 30.0, 60.0] | units.deg,
    longitude_of_ascending_node=[0.001, 0.001, 0.001] | units.deg,
):
    """
    Evolves quadruple system.
    """
    print(longitude_of_ascending_node)

    number_of_bodies = len(masses)
    number_of_binaries = number_of_bodies - 1
    particles, binaries = initialize_multiple_system(
        number_of_bodies,
        masses,
        semimajor_axis,
        eccentricity,
        inclination,
        argument_of_pericenter,
        longitude_of_ascending_node,
    )

    code = Secularmultiple()
    code.particles.add_particles(particles)

    channel_from_particles_to_code = particles.new_channel_to(code.particles)
    channel_from_code_to_particles = code.particles.new_channel_to(particles)
    channel_from_particles_to_code.copy()

    # set up some arrays for plotting
    print_smas_au = [[] for x in range(number_of_binaries)]
    print_rps_au = [[] for x in range(number_of_binaries)]
    print_parent_is_deg = [[] for x in range(number_of_binaries)]
    print_times_myr = []

    time = 0.0 | units.yr
    time_step_output = time_end / float(number_of_output_steps)
    while time <= time_end:
        time += time_step_output
        code.evolve_model(time)

        channel_from_code_to_particles.copy()
        print("=" * 50)
        print(f"time: {time.in_(units.Myr)}")
        print(f"e: {binaries.eccentricity}")
        print(f"i: {binaries.inclination.in_(units.deg)}")
        print(f"AP: {binaries.argument_of_pericenter.in_(units.deg)}")
        print(f"LAN: {binaries.longitude_of_ascending_node.in_(units.deg)}")

        # write to output arrays
        print_times_myr.append(time.value_in(units.Myr))
        for index_binary in range(number_of_binaries):
            print_smas_au[index_binary].append(
                binaries[index_binary].semimajor_axis.value_in(units.au)
            )
            print_rps_au[index_binary].append(
                binaries[index_binary].semimajor_axis.value_in(units.au)
                * (1.0 - binaries[index_binary].eccentricity)
            )
            print_parent_is_deg[index_binary].append(
                np.rad2deg(binaries[index_binary].inclination_relative_to_parent)
            )

    # compute the `canonical' maximum eccentricity/periapsis distance that
    # applies in the quadrupole-order test-particle limit if the `outer' binary
    # is replaced by a point mass
    print(
        inclination[0],
        inclination[2],
        longitude_of_ascending_node[0],
        longitude_of_ascending_node[2],
    )
    i_ac_init = compute_mutual_inclination(
        inclination[0],
        inclination[2],
        longitude_of_ascending_node[0],
        longitude_of_ascending_node[2],
    )
    i_bc_init = compute_mutual_inclination(
        inclination[1],
        inclination[2],
        longitude_of_ascending_node[1],
        longitude_of_ascending_node[2],
    )

    canonical_rp_min_a_au = (
        semimajor_axis[0] * (1.0 - np.sqrt(1.0 - (5.0 / 3.0) * cos(i_ac_init) ** 2))
    ).value_in(units.au)
    canonical_rp_min_b_au = (
        semimajor_axis[1] * (1.0 - np.sqrt(1.0 - (5.0 / 3.0) * cos(i_bc_init) ** 2))
    ).value_in(units.au)

    data = (
        print_times_myr,
        print_smas_au,
        print_rps_au,
        print_parent_is_deg,
        canonical_rp_min_a_au,
        canonical_rp_min_b_au,
    )
    return data


def compute_mutual_inclination(
    inclination_k,
    inclination_l,
    longitude_of_ascending_node_k,
    longitude_of_ascending_node_l,
):
    cos_inclination_rel = cos(inclination_k) * cos(inclination_l) + sin(
        inclination_k
    ) * sin(inclination_l) * cos(
        longitude_of_ascending_node_k - longitude_of_ascending_node_l
    )
    return arccos(cos_inclination_rel)


def plot_function(data):
    (
        print_times_myr,
        print_smas_au,
        print_rps_au,
        print_parent_is_deg,
        canonical_rp_min_a_au,
        canonical_rp_min_b_au,
    ) = data

    number_of_binaries = len(print_smas_au)

    linewidth = 4
    dlinewidth = 2
    fig = plt.figure(figsize=(10, 9))
    plot1 = fig.add_subplot(2, 1, 1, yscale="log")
    plot2 = fig.add_subplot(2, 1, 2)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    labels = ["$A$", "$B$", "$C$"]
    labels_i = ["$i_{AC}$", "$i_{BC}$", None]
    for index_binary in range(number_of_binaries):
        label = labels[index_binary]
        label_i = labels_i[index_binary]
        color = colors[index_binary]

        plot1.plot(
            print_times_myr,
            print_smas_au[index_binary],
            color=color,
            linestyle="dashed",
            linewidth=dlinewidth,
        )
        plot1.plot(
            print_times_myr,
            print_rps_au[index_binary],
            color=color,
            linewidth=linewidth,
            label=label,
        )

        plot2.plot(
            print_times_myr,
            print_parent_is_deg[index_binary],
            color=color,
            linewidth=linewidth,
            label=label_i,
        )

    plot1.axhline(
        y=canonical_rp_min_a_au,
        color=colors[0],
        linestyle="dotted",
        linewidth=dlinewidth,
    )
    plot1.axhline(
        y=canonical_rp_min_b_au,
        color=colors[1],
        linestyle="dotted",
        linewidth=dlinewidth,
    )

    handles, labels = plot1.get_legend_handles_labels()
    plot1.legend(handles, labels, loc="upper right", fontsize=12)

    handles, labels = plot2.get_legend_handles_labels()
    plot2.legend(handles, labels, loc="lower right", fontsize=12)

    plot1.set_xlabel("t [Myr]", fontsize=18)
    plot2.set_xlabel("t [Myr]", fontsize=18)

    plot1.set_ylabel(r"a$_i$ [au]", fontsize=18)
    plot2.set_ylabel(r"i$_{kl}$ [deg]", fontsize=18)

    plot1.set_xlim(0.0, print_times_myr[-1])
    plot2.set_xlim(0.0, print_times_myr[-1])

    plot1.tick_params(axis="both", which="major", labelsize=18)
    plot2.tick_params(axis="both", which="major", labelsize=18)

    fig.savefig("hierarchical_quadruple.pdf")


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=5.0 | units.Myr,
        help="integration time",
    )
    parser.add_argument(
        "-n",
        "--number_of_output_steps",
        type=int,
        default=400,
        help="number of output steps",
    )
    parser.add_argument(
        "--m1",
        "--mass_1",
        type=units.MSun,
        default=1.0 | units.MSun,
        help="mass of object 1",
    )
    parser.add_argument(
        "--m2",
        "--mass_2",
        type=units.MSun,
        default=0.8 | units.MSun,
        help="mass of object 2",
    )
    parser.add_argument(
        "--m3",
        "--mass_3",
        type=units.MSun,
        default=1.1 | units.MSun,
        help="mass of object 3",
    )
    parser.add_argument(
        "--m4",
        "--mass_4",
        type=units.MSun,
        default=0.9 | units.MSun,
        help="mass of object 4",
    )
    parser.add_argument(
        "--aA",
        "--semimajor_axis_a",
        type=units.au,
        default=1.0 | units.au,
        help="semimajor axis of orbit A",
    )
    parser.add_argument(
        "--aB",
        "--semimajor_axis_b",
        type=units.au,
        default=1.2 | units.au,
        help="semimajor axis of orbit B",
    )
    parser.add_argument(
        "--aC",
        "--semimajor_axis_c",
        type=units.au,
        default=100.0 | units.au,
        help="semimajor axis of orbit C (the 'superorbit')",
    )
    parser.add_argument(
        "--eA",
        "--eccentricity_a",
        type=float,
        default=0.1,
        help="eccentricity of orbit A",
    )
    parser.add_argument(
        "--eB", 
        "--eccentricity_b",
        type=float,
        default=0.1,
        help="eccentricity of orbit B",
    )
    parser.add_argument(
        "--eC",
        "--eccentricity_c",
        type=float,
        default=0.3,
        help="eccentricity of orbit C (the 'superorbit')",
    )
    parser.add_argument(
        "--iA",
        "--inclination_a",
        type=float,
        default=75.0 | units.deg,
        help="inclination of orbit A in degrees",
    )
    parser.add_argument(
        "--iB",
        "--inclination_b",
        type=float,
        default=80.0 | units.deg,
        help="inclination of orbit B in degrees",
    )
    parser.add_argument(
        "--iC",
        "--inclination_c",
        type=float,
        default=0.001 | units.deg,
        help="inclination of orbit C (the 'superorbit') in degrees",
    )
    parser.add_argument(
        "--ApA",
        "--argument_of_periapsis_a",
        type=float,
        default=10.0 | units.deg,
        help="argument of periapsis of orbit A in degrees",
    )
    parser.add_argument(
        "--ApB",
        "--argument_of_periapsis_b",
        type=float,
        default=30.0 | units.deg,
        help="argument of periapsis of orbit B in degrees",
    )
    parser.add_argument(
        "--ApC",
        "--argument_of_periapsis_c",
        type=float,
        default=60.0 | units.deg,
        help=("argument of periapsis of orbit C (the 'superorbit') in degrees"),
    )
    parser.add_argument(
        "--LANA",
        "--longitude_of_ascending_node_a",
        type=float,
        default=0.001 | units.deg,
        help=("longitude of the ascending node of orbit A in degrees"),
    )
    parser.add_argument(
        "--LANB",
        "--longitude_of_ascending_node_b",
        type=float,
        default=0.001 | units.deg,
        help=("longitude of the ascending node of orbit B in degrees"),
    )
    parser.add_argument(
        "--LANC",
        "--longitude_of_ascending_node_c",
        type=float,
        default=0.001 | units.deg,
        help=(
            "longitude of the ascending node of orbit C (the 'superorbit') "
            "in degrees"
        ),
    )

    return parser


def hierarchical_quadruple(
    number_of_output_steps=400,
    time_end=5.0 | units.Myr,
    masses=[1.0, 0.8, 1.1, 0.9] | units.MSun,
    semimajor_axis=[1.0, 1.2, 100.0] | units.au,
    eccentricity=[0.1, 0.1, 0.3],
    inclination=[75.0, 80.0, 0.001] | units.deg,
    argument_of_pericenter=[10.0, 30.0, 60.0] | units.deg,
    longitude_of_ascending_node=[0.001, 0.001, 0.001] | units.deg,
):
    data = evolve_quadruple(
        number_of_output_steps=number_of_output_steps,
        time_end=time_end,
        masses=masses,
        semimajor_axis=semimajor_axis,
        eccentricity=eccentricity,
        inclination=inclination,
        argument_of_pericenter=argument_of_pericenter,
        longitude_of_ascending_node=longitude_of_ascending_node,
    )
    plot_function(data)


def main(**kwargs):
    hierarchical_quadruple(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    masses = as_vector_quantity(
        [
            arguments.mass_1,
            arguments.mass_2,
            arguments.mass_3,
            arguments.mass_4,
        ]
    )
    semi_major_axes = as_vector_quantity(
        [
            arguments.semimajor_axis_a,
            arguments.semimajor_axis_b,
            arguments.semimajor_axis_c,
        ]
    )
    eccentricities = as_vector_quantity(
        [
            arguments.eccentricity_a,
            arguments.eccentricity_b,
            arguments.eccentricity_c,
        ]
    )
    inclinations = as_vector_quantity(
        [
            arguments.inclination_a,
            arguments.inclination_b,
            arguments.inclination_c,
        ]
    )
    arguments_of_periapsis = as_vector_quantity(
        [
            arguments.argument_of_periapsis_a,
            arguments.argument_of_periapsis_b,
            arguments.argument_of_periapsis_c,
        ]
    )
    longitudes_of_ascending_node = as_vector_quantity(
        [
            arguments.longitude_of_ascending_node_a,
            arguments.longitude_of_ascending_node_b,
            arguments.longitude_of_ascending_node_c,
        ]
    )
    kwargs = dict(
        masses=masses,
        semimajor_axis=semi_major_axes,
        eccentricity=eccentricities,
        inclination=inclinations,
        argument_of_pericenter=arguments_of_periapsis,
        longitude_of_ascending_node=longitudes_of_ascending_node,
    )
    main(**kwargs)
