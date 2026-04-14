import sys
import numpy as np

import matplotlib.pyplot as plt

from amuse.units import units, nbody_system
from amuse.units.optparse import OptionParser
from amuse.datamodel import Particle
from amuse.community.mesa import Mesa
from amuse.community.gadget2 import Gadget2
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.relax_sph import relax


def evolve_star(mass, age):
    stellar = Mesa()
    stellar.parameters.metallicity = 0.02
    star = stellar.particles.add_particle(Particle(mass=mass))
    stellar.evolve_model(age)
    return stellar


def star_to_sph(mass, age, omega, number_of_sph_particles, core_mass):
    stellar = evolve_star(mass, age)
    # new_age = stellar.model_time
    star = stellar.particles[0]
    if core_mass > 0 | units.MSun:
        with_core_particle = True
        if star.core_mass > 0 | units.MSun:
            target_core_mass = star.core_mass
        else:
            target_core_mass = core_mass
    else:
        with_core_particle = False
        target_core_mass = 0 | units.MSun

    m_sph = (star.mass-target_core_mass)/float(number_of_sph_particles)
    print(f"N_sph= {number_of_sph_particles} {m_sph.in_(units.MSun)}")
    print(f"Target core mass: {target_core_mass}")
    model = convert_stellar_model_to_SPH(
        star,
        number_of_sph_particles,
        with_core_particle=with_core_particle,
        target_core_mass=target_core_mass,
        do_store_composition=False,
        base_grid_options={"type": "fcc"},
    )
    print(f"Final star: {star}")
    age = stellar.model_time
    stellar.stop()

    core = model.core_particle
    print(f"core {core}")
    if core is None:
        print("Make zero mass core")
        core = Particle(mass=0 | units.MSun, radius=1 | units.RSun)
        core.position = (0, 0, 0) | units.AU
        core.velocity = (0, 0, 0) | units.kms
    sph_star = model.gas_particles
    if omega != 0 | units.day**-1:
        print(f"Add Spin to star: {omega}")
        sph_star.add_spin(omega)

    relaxed_sph_star, core_particles = relax_sph_realization(sph_star)
    return relaxed_sph_star, core_particles, age


def relax_sph_realization(sph_star):
    dynamical_timescale = sph_star.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1 | units.RSun)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.parameters.time_max = 3 * dynamical_timescale
    hydro.parameters.max_size_timestep = dynamical_timescale / 100
    hydro.parameters.time_limit_cpu = 1.0 | units.Gyr

    t_end_in_t_dyn = 2.5  # Relax for this many dynamical timescales
    t_end = t_end_in_t_dyn * sph_star.dynamical_timescale(mass_fraction=0.9)
    n_steps = 250
    velocity_damp_factor = (
        1.0 - (2.0*np.pi*t_end_in_t_dyn)/n_steps  # Critical damping
    )

    hydro.gas_particles.add_particles(sph_star)
    to_framework = hydro.gas_particles.new_channel_to(sph_star)

    # stellar_to_framework = stellar.particles.new_channel_to(stars)
    # if core.mass>0|units.MSun:
    #    hydro.dm_particles.add_particles(core)
    for i_step, time in enumerate(t_end * np.linspace(1.0/n_steps, 1.0, n_steps)):
        hydro.evolve_model(time)
        to_framework.copy()
        hydro.gas_particles.velocity = velocity_damp_factor * hydro.gas_particles.velocity

    return hydro.gas_particles, hydro.dm_particles


def get_density_profile(
    code=Mesa, mass=1.0 | units.MSun, metallicity=0.02, model_time=1 | units.yr
):
    stellar = code()
    stellar.parameters.metallicity = metallicity
    stellar.particles.add_particle(Particle(mass=mass))
    stellar.evolve_model(model_time)
    print("Nzones=", stellar.particles.get_number_of_zones())
    radius = stellar.particles[0].get_radius_profile()
    rho = stellar.particles[0].get_density_profile()
    stellar.stop()
    return radius, rho


def plot_star_in_sph(star, core, mass):
    color = plt.rcParams['axes.prop_cycle'].by_key()['color']

    x_label = r"$R$ [R$_\odot$]"
    y_label = "$\\rho$ [g/cm$^{3}$]"
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_yscale('log')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    plt.xlim(0, 2)
    plt.ylim(1.e-9, 1.e+2)

    radius = np.array([])
    density = np.array([])
    for si in star:
        radius = np.append(radius, si.position.length().value_in(units.RSun))
        density = np.append(density, si.density.value_in(units.g/units.cm**3))
    radius, density = zip(*sorted(zip(radius, density)))
    number_of_particles = len(star)
    filename = (
        f"fig_star_in_sph_N{number_of_particles}M{(mass.value_in(units.MSun))}MSun"
    )
    with open(f"{filename}.npy", 'wb') as save_file:
        np.save(save_file, radius)
        np.save(save_file, density)

    print(f"r= {len(radius)} {len(density)}")
    plt.scatter(radius, density, c=color[0], s=2)

    metallicity = 0.02
    radius, density = get_density_profile(Mesa, mass, metallicity)
    plt.plot(
        radius.value_in(units.RSun),
        density.value_in(units.g/units.cm**3),
        label="MESA", c=color[1],
    )
    plt.legend(loc="lower right")
    plt.semilogy()
    plt.savefig(filename+'.pdf')
    plt.show()


def new_option_parser():
    result = OptionParser()
    result.add_option("-N",
                      dest="number_of_sph_particles", type="int",
                      default=1000,
                      help="number of SPH particles[%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="mass", type="float",
                      default=1.0 | units.MSun,
                      help="stellar mass [%default]")
    result.add_option("--mcore", unit=units.MSun,
                      dest="core_mass", type="float",
                      default=-1.0 | units.MSun,
                      help="core mass [%default]")
    result.add_option("-t", unit=units.Myr,
                      dest="age", type="float",
                      default=1 | units.yr,
                      help="stellar age [%default]")
    result.add_option("-o", unit=units.day**-1,
                      dest="omega", type="float",
                      # default=24.47 | units.day**-1,
                      default=0 | units.day**-1,
                      help="stellar rotation [%default]")
    return result


def main():
    o, arguments = new_option_parser().parse_args()
    if o.mass <= o.core_mass:
        print("Core mass should not exceed the stellar mass")
        sys.exit()
    star, core, age = star_to_sph(o.mass, o.age, o.omega, o.number_of_sph_particles, o.core_mass)
    print(f"age={age.in_(units.Myr)}")
    print(star)
    print(star.mass.sum().in_(units.MSun), core.mass)
    filename = f'Hydro_M{o.mass.value_in(units.MSun):2.2d}MSun.amuse'
    print(filename)

    plot_star_in_sph(star, core, o.mass)
#    write_set_to_file(star, filename, format='hdf5', append_to_file=False)
#    write_set_to_file(core, filename, format='hdf5', append_to_file=True)


if __name__ == "__main__":
    main()
