from amuse.datamodel import Particles
from amuse.units import units

from amuse.community.mesa import Mesa
from amuse.community.mmams import Mmams

from amuse.couple.collision_handler import CollisionHandler

from amuse.support.console import set_printing_strategy

from matplotlib import pyplot as plt

DEFAULT_OPTIONS = {
    "redirection": "none",
}
MMAMS_OPTIONS = {**DEFAULT_OPTIONS}
MESA_OPTIONS = {
    **DEFAULT_OPTIONS,
    "version": "2208",
}

def evolve_single_star(mass, time_end):
    star = Particles(1)
    star.mass = mass
    stellar_evolution = Mesa(**MESA_OPTIONS)
    stellar_evolution.parameters.metallicity = 0.0142

    stellar_evolution.particles.add_particles(star)
    time = [] | units.Myr
    stellar_type = []
    mass = [] | units.MSun
    radius = [] | units.RSun
    temperature = [] | units.K
    luminosity = [] | units.LSun
    while stellar_evolution.model_time < time_end:
        stellar_evolution.evolve_model()
        time.append(stellar_evolution.model_time)
        stellar_type.append(stellar_evolution.particles[0].stellar_type)
        mass.append(stellar_evolution.particles[0].mass)
        radius.append(stellar_evolution.particles[0].radius)
        temperature.append(stellar_evolution.particles[0].temperature)
        luminosity.append(stellar_evolution.particles[0].luminosity)
        print(
            f"Time={time[-1]} {stellar_type[-1]} {mass[-1]} {radius[-1]} "
            f"{temperature[-1].in_(units.K)} {luminosity[-1].in_(units.LSun)}"
        )
        if stellar_type[-1] >= 4 | units.stellar_type:
            break

    stellar_evolution.stop()
    return time, stellar_type, mass, radius, temperature, luminosity

def print_stars(stellar_evolution):
    print(
        "Primary:   Time=", stellar_evolution.model_time,
        stellar_evolution.particles[0].mass,
        stellar_evolution.particles[0].radius,
        stellar_evolution.particles[0].temperature,
        stellar_evolution.particles[0].luminosity
    )
    print(
        "Secondary: Time=", stellar_evolution.model_time,
        stellar_evolution.particles[1].mass,
        stellar_evolution.particles[1].radius,
        stellar_evolution.particles[1].temperature,
        stellar_evolution.particles[1].luminosity
    )

# #BOOKLISTSTART1# #
def merge_two_stars_and_evolve(
    mass_primary, mass_secondary, time_collision, time_end
):
    stars = Particles(2)
    stars.mass = [
        mass_primary.value_in(units.MSun),
        mass_secondary.value_in(units.MSun)
    ] | units.MSun
    stellar_evolution = Mesa(**MESA_OPTIONS)
    stellar_evolution.parameters.metallicity = 0.0142
    stellar_evolution.particles.add_particles(stars)
    while stellar_evolution.model_time < time_collision:
        stellar_evolution.evolve_model()
        print_stars(stellar_evolution)
    n_shell = min(
        stellar_evolution.particles[0].get_number_of_zones(),
        stellar_evolution.particles[1].get_number_of_zones(),
    )
    merger_code = Mmams(**MMAMS_OPTIONS)
    merger_code.parameters.target_n_shells = n_shell
    merger_code.parameters.dump_mixed_flag = True
    merger_code.parameters.do_shock_heating_flag = True
    merger_code.commit_parameters()
# #BOOKLISTSTOP1# #

# #BOOKLISTSTART2# #
    handler = CollisionHandler(
        merger_code,
        stellar_evolution_code=stellar_evolution,
    )
    merger_product = handler.handle_collision(
        stellar_evolution.particles[0],
        stellar_evolution.particles[1],
    )
    merged = stellar_evolution.particles[0]
# #BOOKLISTSTOP2# #

    print("Stars merged:", merged)
    time = [] | units.Myr
    stellar_type = []
    mass = [] | units.MSun
    radius = [] | units.RSun
    temperature = [] | units.K
    luminosity = [] | units.LSun
    
# #BOOKLISTSTART3# #
    stellar_evolution.evolve_model(keep_synchronous=True)
    p = stellar_evolution.particles[0]
    while stellar_evolution.model_time < time_end:
        stellar_evolution.evolve_model()
# #BOOKLISTSTOP3# #

        time.append(stellar_evolution.model_time)
        stellar_type.append(p.stellar_type)
        mass.append(p.mass)
        radius.append(p.radius)
        temperature.append(p.temperature)
        luminosity.append(p.luminosity)

# #BOOKLISTSTART4# #
        print(
            f"Time={stellar_evolution.model_time} {p.stellar_type} "
            f"{p.mass} {p.radius} {p.temperature} {p.luminosity}"
        )
        if p.stellar_type >= 4 | units.stellar_type:
            break
    merger_code.stop()
    stellar_evolution.stop()
# #BOOKLISTSTOP4# #

    return time, stellar_type, mass, radius, temperature, luminosity

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--tcoll", unit=units.Myr,
                      dest="tcoll", type="float",
                      default=150 | units.Myr,
                      help="moment of collision [%default]")
    result.add_option("--tend", unit=units.Myr,
                      dest="tend", type="float",
                      default=2 | units.Gyr,
                      help="evolution after the collision [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mprim", type="float",
                      default=3 | units.MSun,
                      help="Primary ZAMS mass [%default]")
    result.add_option("-m", unit=units.MSun,
                      dest="Msec", type="float",
                      default=1 | units.MSun,
                      help="Secondary ZAMS mass [%default]")

    return result


def main():
    # High-level structure of merge_two_stars_and_evolve.py and
    # merge_two_stars_sph_evolve.py are designed to be identical.

    set_printing_strategy(
        "custom",  # nbody_converter=converter,
        precision=11, prefix="",
        separator=" [",
        suffix="]",
        preferred_units=(
            units.MSun, units.K, units.LSun, units.RSun,
            units.Myr,
        )
    )

    o, arguments = new_option_parser().parse_args()
    mass_primary = o.Mprim
    mass_secondary = o.Msec
    time_end = o.tend
    time_collision = o.tcoll

    color = plt.rcParams['axes.prop_cycle'].by_key()['color']

    x_label = "T [K]"
    y_label = r"L [$L_\odot$]"
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    print(f"Evolve single star of mass {mass_primary.in_(units.MSun)}")
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(mass_primary, time_end)
    plt.plot(
        temperature.value_in(units.K),
        luminosity.value_in(units.LSun),
        c=color[1], lw=2, zorder=1,
        )
    plt.scatter(
        temperature[0].value_in(units.K),
        luminosity[0].value_in(units.LSun),
        c=color[1], s=150, marker='^',
        edgecolor='k', zorder=2,
    )

    time_ms = 0 | units.Myr
    for i, stellar_type in enumerate(stp):
        if stellar_type < 2 | units.stellar_type:
            time_ms = time[i]
    if time_ms <= 1 | units.Myr:
        time_ms = 10 | units.Myr
    print(f"Main-sequence lifetime = {time_ms.in_(units.Myr)}")

    tcoll = 0.5*time_ms
    icoll = 0
    for i in range(len(stp)):
        if time[i] <= tcoll:
            icoll = i
    plt.scatter(
        temperature[icoll].value_in(units.K),
        luminosity[icoll].value_in(units.LSun),
        c=color[2], s=150, marker='o',
        edgecolor='k', zorder=2,
    )

    print(
        f"Evolve single star of mass {(mass_primary+mass_secondary)}"
    )
    time, stp, mass, radius, temperature, luminosity \
        = evolve_single_star(mass_primary+mass_secondary, time_end)
    plt.plot(
        temperature.value_in(units.K),
        luminosity.value_in(units.LSun),
        c=color[0], lw=2, zorder=1,
    )
    plt.scatter(
        temperature[0].value_in(units.K),
        luminosity[0].value_in(units.LSun),
        c=color[0], s=150, marker='^',
        edgecolor='k', zorder=2,
    )

    print(f"Evolve two single stars and collide at {tcoll.in_(units.Myr)}")
    time, stp, mass, radius, temperature, luminosity \
        = merge_two_stars_and_evolve(
            mass_primary, mass_secondary, time_collision, time_end
        )
    plt.plot(
        temperature.value_in(units.K),
        luminosity.value_in(units.LSun),
        c=color[2], ls="--", lw=3, zorder=1,
    )
    plt.scatter(
        temperature[0].value_in(units.K),
        luminosity[0].value_in(units.LSun),
        c=color[2], s=150, marker='o',
        edgecolor='k', zorder=3,
    )

    mass_merger = mass[0]
    print(f"Evolve single star of mass {mass_merger}")
    time, stp, mass, radius, temperature, luminosity = \
        evolve_single_star(mass_merger, time_end)
    plt.plot(
        temperature.value_in(units.K),
        luminosity.value_in(units.LSun),
        c=color[3], lw=2, zorder=1,
        )
    plt.scatter(
        temperature[0].value_in(units.K),
        luminosity[0].value_in(units.LSun),
        c=color[3], s=150, marker='^',
        edgecolor='k', zorder=2,
    )

    ax = plt.gca()
    ax.tick_params(axis='both', which='both', direction='in')
    ax.invert_xaxis()

    save_file = 'merge_two_stars_and_evolve.pdf'
    plt.savefig(save_file)
    print(f'Saved figure in file {save_file}')
    plt.show()


if __name__ in ('__main__', '__plot__'):
    main()
