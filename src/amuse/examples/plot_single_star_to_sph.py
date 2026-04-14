import sys
if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files
import numpy as np
from matplotlib import pyplot as plt

from .single_star_to_sph import get_density_profile


def plot_single_star_to_sph():
    datadir = files("amuse.examples.data")
    with open(datadir.joinpath("fig_star_in_sph_N10000M1.0MSun.npy"), "rb") as f:
        r1MSun = np.load(f)
        rho1MSun = np.load(f)

    with open(datadir.joinpath("fig_star_in_sph_N10000M3.0MSun.npy"), "rb") as f:
        r3MSun = np.load(f)
        rho3MSun = np.load(f)

    with open(datadir.joinpath("fig_star_in_sph_N10000M10.0MSun.npy"), "rb") as f:
        r10MSun = np.load(f)
        rho10MSun = np.load(f)
    
    color = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    x_label = r"$R$ [R$_\odot$]"
    y_label = r"$\\rho$ [g/cm$^{3}$]"
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_yscale('log')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    plt.xlim(1.e-2, 5)
    plt.ylim(2.e-3, 1.e+2)

    plt.scatter(r1MSun, rho1MSun, c=color[0], s=2)
    plt.scatter(r3MSun, rho3MSun, c=color[1], s=2)
    plt.scatter(r10MSun, rho10MSun, c=color[2], s=2)

    z = 0.02
    mass = 1|units.MSun
    r, rho = get_density_profile(MESA, mass, z)
    plt.plot(r.value_in(units.RSun),
             rho.value_in(units.g/units.cm**3),
             label="$1 M_\odot$ ZAMS", c=color[0])
    r, rho = get_density_profile(MESA, 3|units.MSun, z)
    plt.plot(r.value_in(units.RSun),
             rho.value_in(units.g/units.cm**3),
             label="$3 M_\odot$", c=color[1])
    r, rho = get_density_profile(MESA, 10|units.MSun, z)
    plt.plot(r.value_in(units.RSun),
             rho.value_in(units.g/units.cm**3),
             label="$10 M_\odot$", c=color[2])
    plt.legend(loc="lower left")
    #plt.semilogy()
    plt.loglog()
    filename = "fig_star_in_sph_N1e4at0Myr"
    plt.savefig(f"{filename}.pdf")


def main():
    plot_single_star_to_sph()


if __name__ == "__main__":
    main()
