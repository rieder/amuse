import argparse
from amuse.datamodel import Particles
from amuse.units import units, nbody_system
from amuse.ic.gasplummer import new_plummer_gas_model


def accrete_particles(sinks, gas):
    removed_particles = Particles()
    for s in sinks:
        xs, ys, zs = s.x, s.y, s.z
        radius_squared = s.radius**2
        insink = gas.select_array(
            lambda x, y, z: (x - xs) ** 2 + (y - ys) ** 2 + (z - zs) ** 2
            < radius_squared,
            ["x", "y", "z"],
        )
        if len(insink) == 0:
            return insink

        cm = s.position * s.mass
        p = s.velocity * s.mass
        s.mass += insink.total_mass()
        s.position = (cm + insink.center_of_mass() * insink.total_mass()) / s.mass
        s.velocity = (p + insink.total_momentum()) / s.mass
        removed_particles.add_particles(insink)
    return removed_particles


def hydro_sink_particles(
    number_of_particles=100,
    mass_cloud=1 | units.MSun,
    radius_cloud=100 | units.au,
    radius_sink=1 | units.au,
):
    converter = nbody_system.nbody_to_si(mass_cloud, radius_cloud)
    bodies = new_plummer_gas_model(number_of_particles, convert_nbody=converter)

    sink = Particles(1)
    sink.mass = 0 | units.MSun
    sink.radius = radius_sink
    sink.position = (0, 0, 0) | units.au
    sink.velocity = (0, 0, 0) | units.kms

    accreted = accrete_particles(sink, bodies)
    print(f"Particles in sink: N={len(accreted)} M={sink.mass}")
    print(f"sink position={sink.position.as_quantity_in(units.au)}")
    print(f"sink velocity={sink.velocity.as_quantity_in(units.kms)}")


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-N",
        "--number_of_particles",
        type=int,
        default=100,
        help="number of sph particles",
    )
    parser.add_argument(
        "-M",
        "--mass_cloud",
        type=units.MSun,
        default=1 | units.MSun,
        help="Mass of molecular cloud",
    )
    parser.add_argument(
        "-R",
        "--radius_cloud",
        type=units.au,
        default=100 | units.au,
        help="Virial radius of cloud",
    )
    parser.add_argument(
        "-r",
        "--radius_sink",
        type=units.au,
        default=100 | units.au,
        help="Radius of the sink",
    )
    return parser


def main(**kwargs):
    hydro_sink_particles(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
