"""
Creates a helium burning star from the inner shells of a main sequence star.
"""

import matplotlib.pyplot as plt

from amuse.units import units
from amuse.datamodel import Particle
from amuse.community.mesa import Mesa
from amuse.plot import loglog, xlabel, ylabel, scatter


def make_helium_star():
    temperatures_original = [] | units.K
    luminosities_original = [] | units.LSun
    temperatures_helium = [] | units.K
    luminosities_helium = [] | units.LSun

    star = Particle()
    star.mass = 12.0 | units.MSun
    stop_radius = 90 | units.RSun  # Avoid the model going through a blue loop

    stellar_evolution = Mesa()
    se_star = stellar_evolution.particles.add_particle(star)

    # This method is only really useful if you need a helium star from an
    # existing model.
    # With the new interface you can do:
    # stellar_evolution.pure_he_stars.add_particles(star)
    # (which is much faster).
    # TODO check if this is still true for Mesa 15140 -- SR

    print(
        f"Evolving a {star.mass} star with {stellar_evolution.__class__.__name__}"
        f" until its radius exceeds {stop_radius}"
    )
    while se_star.radius < stop_radius:
        se_star.evolve_one_step()
        temperatures_original.append(se_star.temperature)
        luminosities_original.append(se_star.luminosity)

    # #BOOKLISTSTART1# #
    number_of_zones = se_star.get_number_of_zones()
    composition = se_star.get_chemical_abundance_profiles()
    # first zone with X > 1.0e-9
    index = (composition[0] > 1.0e-9).nonzero()[0][0]

    print(
        f"Creating helium star, from the inner {index} (out of {number_of_zones})"
        f"shells."
    )
    helium_star_in_code = stellar_evolution.new_particle_from_model(
        {
            "mass": (se_star.get_cumulative_mass_profile() * se_star.mass)[:index],
            "radius": se_star.get_radius_profile()[:index],
            "rho": se_star.get_density_profile()[:index],
            "temperature": se_star.get_temperature_profile()[:index],
            "luminosity": se_star.get_luminosity_profile()[:index],
            "X_H": composition[0][:index],
            "X_He": composition[1][:index] + composition[2][:index],
            "X_C": composition[3][:index],
            "X_N": composition[4][:index],
            "X_O": composition[5][:index],
            "X_Ne": composition[6][:index],
            "X_Mg": composition[7][:index],
            "X_Si": composition[7][:index] * 0.0,
            "X_Fe": composition[7][:index] * 0.0,
        },
        0.0 | units.Myr,
    )
    # #BOOKLISTSTOP1# #
    helium_star = Particle()

    print(
        "\nStar properties before helium star evolution:\n", stellar_evolution.particles
    )
    for i in range(1000):
        helium_star_in_code.evolve_one_step()
        temperatures_helium.append(helium_star_in_code.temperature)
        luminosities_helium.append(helium_star_in_code.luminosity)
    print(
        "\nStar properties after helium star evolution:\n", stellar_evolution.particles
    )

    stellar_evolution.stop()
    # FIXME: save this in an AMUSE particle instead
    return (
        temperatures_original,
        luminosities_original,
        temperatures_helium,
        luminosities_helium,
    )


def plot_tracks(
    temperatures_original,
    luminosities_original,
    temperatures_helium,
    luminosities_helium,
):
    """
    Plots tracks for temperature and luminosity of the original star and the
    helium star.
    """

    x_label = "T [K]"
    y_label = r"L [$L_\odot$]"
    figure = plt.figure()
    ax = figure.add_subplot(1, 1, 1)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # single_frame(x_label, y_label, logx=True, logy=True, xsize=14, ysize=10)
    # colors = get_distinct(2)
    prop_cycle = plt.rcParams["axes.prop_cycle"]
    colors = prop_cycle.by_key()["color"]

    loglog(
        temperatures_original, luminosities_original, label="progenitor", c=colors[0]
    )
    loglog(temperatures_helium, luminosities_helium, label="helium star", c=colors[1])
    scatter(
        temperatures_helium[-1], luminosities_helium[-1], marker="*", s=400, c=colors[1]
    )
    xlabel("Effective Temperature")
    ylabel("Luminosity")
    plt.xlim(1.0e5, 4000)
    plt.ylim(1.0, 1.0e5)
    plt.legend(loc=4, fontsize=24)

    save_file = "HertzsprungRussell_HeliumStar.png"
    plt.savefig(save_file)
    print(f"Saved figure in file {save_file}\n")
    plt.show()


def main():
    (
        temperatures_original,
        luminosities_original,
        temperatures_helium,
        luminosities_helium,
    ) = make_helium_star()
    plot_tracks(
        temperatures_original,
        luminosities_original,
        temperatures_helium,
        luminosities_helium,
    )


if __name__ == "__main__":
    main()
