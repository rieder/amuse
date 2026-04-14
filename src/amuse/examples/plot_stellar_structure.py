import matplotlib.pyplot as plt
import seaborn as sns

from amuse.datamodel import Particle
from amuse.units import units
from amuse.community.mesa import Mesa


def stellar_density_profile_at_time(mass, time):
    stellar = Mesa()
    stellar.parameters.metallicity = 0.02
    star = stellar.particles.add_particle(Particle(mass=mass))

    stellar.evolve_model(time)

    density_profile = star.get_density_profile(star.get_number_of_zones())
    luminosity_profile = star.get_luminosity_profile(
        star.get_number_of_zones()
    )
    mass_profile = star.get_cumulative_mass_profile(
        star.get_number_of_zones()
    ) * mass
    stellar.stop()

    return mass_profile, luminosity_profile, density_profile


def plot_stellar_structure():
    figure = plt.figure()
    ax = figure.add_subplot(111)

    ax.set_xlabel(r"L [$L_\odot$]")
    ax.set_ylabel("density [$g/cm^{3}$]")
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1.e-9, 1.e+5)
    ax.set_ylim(1.e-10, 1e+4)

    masses = [1, 1, 5, 5, 10, 10] | units.MSun
    age = [1, 10000, 1, 90, 1, 20] | units.Myr
    colors = [
        sns.color_palette()[0],
        sns.color_palette()[1],
        sns.color_palette()[2],
    ]
    color = []
    line_styles = []
    for color_i in colors:
        color.append(color_i)
        color.append(color_i)
        line_styles.append("-")
        line_styles.append("--")
    symbol = ['v', 'v', 'o', 'o', '^', '^']
    for i, mass in enumerate(masses):
        dm = 0.2*mass
        time = age[i]
        mass_profile, luminosity_profile, density_profile = \
            stellar_density_profile_at_time(mass, time)
        plt.plot(
            luminosity_profile.value_in(units.LSun),
            density_profile.value_in(units.g/units.cm**3),
            lw=4, color=color[i], label=r'$10M_\odot$', ls=line_styles[i])
        mlim = dm
        for j, mass_shell in enumerate(mass_profile):
            if mass_shell > mlim:
                mlim += dm
                print(j, len(luminosity_profile), len(density_profile))
                plt.scatter(
                    luminosity_profile[j].value_in(units.LSun),
                    density_profile[j].value_in(units.g/units.cm**3),
                    color=color[i], s=100, marker=symbol[i],
                    lw=0
                )
        plt.scatter(
            luminosity_profile[-1].value_in(units.LSun),
            density_profile[-1].value_in(units.g/units.cm**3),
            color=color[i], s=100, marker=symbol[i], lw=0)

    save_file = 'fig_1_5_10_MSun_stellar_core_luminosity.pdf'
    plt.savefig(save_file)
    print(f'\nSaved figure in file {save_file}\n')
    # plt.show()


if __name__ == "__main__":
    plot_stellar_structure()
