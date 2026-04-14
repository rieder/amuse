import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particle, Particles
from amuse.community.fi import Fi
from amuse.community.sphray import Sphray
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.plot import sph_particles_plot


def plot_projection(source, ism, temperature):
    x = ism.x.value_in(units.pc)
    y = ism.y.value_in(units.pc)
    z = ism.z.value_in(units.pc)
    T = temperature.value_in(units.K)

    fig, ax = plt.subplots(figsize=(6, 6))

    cmb = ax.scatter(x, y, c=T, alpha=0.04, s=1000, cmap="rainbow", lw=0)
    ax.scatter(
        source.x.value_in(units.pc),
        source.y.value_in(units.pc),
        c="k",
        s=100,
        marker="*",
    )
    ax.set_aspect(1.0)
    ax.set_xlabel("x [pc]")
    ax.set_ylabel("y [pc]")
    divider = make_axes_locatable(ax)
    ax_histx = divider.append_axes("top", 1.2, pad=0.0, sharex=ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.0, sharey=ax)

    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histx.set_ylabel("z [pc]")
    ax_histy.set_xlabel("z [pc]")
    ax_histy.yaxis.set_tick_params(labelleft=False)

    ax_histx.scatter(x, z, c=T, alpha=0.04, s=100, cmap="rainbow", lw=0)
    ax_histx.scatter(
        source.x.value_in(units.pc),
        source.z.value_in(units.pc),
        c="k",
        s=100,
        marker="*",
    )
    ax_histy.scatter(z, y, c=T, alpha=0.04, s=100, cmap="rainbow", lw=0)
    ax_histy.scatter(
        source.z.value_in(units.pc),
        source.y.value_in(units.pc),
        c="k",
        s=100,
        marker="*",
    )
    ax_histx.set_yticks([-2, 0, 2])
    ax_histy.set_xticks([-2, 0, 2])
    ax.set_xticks([-4, -2, 0, 2, 4])
    ax.set_yticks([-4, -2, 0, 2, 4])

    ax.set_xlim((-5, 5))
    ax.set_xlim((-5, 5))
    ax_histx.set_xlim((-5, 5))
    ax_histx.set_ylim((-3, 3))
    ax_histy.set_xlim((-3, 3))
    ax_histy.set_ylim((-5, 5))

    plt.savefig("fig_projected_disk_and_cloud.pdf")


def mu(X=None, Y=0.25, Z=0.02, x_ion=0.1):
    """
    Compute the mean molecular weight in kg (the average weight of
    particles in a gas) X, Y, and Z are the mass fractions of
    Hydrogen, of Helium, and of metals, respectively.  x_ion is the
    ionisation fraction (0 < x_ion < 1), 1 means fully ionised.
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception(
            "Error in calculating mu: mass " + "fractions do not sum to 1.0"
        )
    return constants.proton_mass / (
        X * (1.0 + x_ion) + Y * (1.0 + 2.0 * x_ion) / 4.0 + Z * x_ion / 2.0
    )


def figure_frame(x_label, y_label, xsize=6, ysize=5):
    figure = plt.figure(figsize=(xsize, ysize))
    ax = figure.add_subplot(1, 1, 1)
    ax.minorticks_on()  # switch on the minor ticks
    ax.locator_params(nbins=3)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    return figure, ax


def plot_gas_temperature(radiative):
    source = radiative.src_particles
    ism = radiative.gas_particles

    mux = mu()
    T = mux / constants.kB * ism.u

    figure = plt.figure(figsize=(10, 10))
    plt.rcParams.update({"font.size": 18})

    ax = figure.add_subplot(projection="3d")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    plt.ticklabel_format(axis="both", style="sci", scilimits=(4, 4))
    ax.scatter(
        ism.x.value_in(units.pc),
        ism.y.value_in(units.pc),
        ism.z.value_in(units.pc),
        s=1000,
        c=T.value_in(units.K),
        alpha=0.04,
        cmap="rainbow",
        lw=0,
    )
    ax.set_facecolor("white")
    plt.savefig("rad_minimal_01.pdf")

    plot_projection(source, ism, T)

    figure = plt.figure(figsize=(10, 10))
    plt.rcParams.update({"font.size": 18})
    ax = plt.gca()
    ax.minorticks_on()  # switch on the minor ticks
    ax.locator_params(nbins=3)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.scatter(
        ism.x.value_in(units.pc),
        ism.y.value_in(units.pc),
        s=1000,
        c=T.value_in(units.K),
        alpha=0.04,
    )
    ax.set_facecolor("white")
    plt.savefig("rad_minimal_02.pdf")

    figure = plt.figure(figsize=(10, 10))
    ax = plt.gca()
    ax.minorticks_on()  # switch on the minor ticks
    ax.locator_params(nbins=3)
    ax.set_xlabel("X")
    ax.set_ylabel("Z")
    plt.ticklabel_format(axis="both", style="sci", scilimits=(4, 4))
    ax.scatter(
        ism.x.value_in(units.pc),
        ism.z.value_in(units.pc),
        s=1000,
        c=T.value_in(units.K),
        alpha=0.04,
    )
    ax.set_facecolor("white")
    plt.savefig("rad_minimal_03.pdf")

    fig, ax = figure_frame("X", "Y", xsize=12, ysize=10)
    cbar = plt.scatter(
        ism.x.value_in(units.pc),
        ism.y.value_in(units.pc),
        ism.z.value_in(units.pc),
        c=T.value_in(units.K),
        # s=1000,
        alpha=0.04,
    )
    plt.colorbar(cbar)
    plt.scatter(source.x.value_in(units.pc), source.y.value_in(units.pc), s=100, c="r")
    plt.savefig("rad_minimal_04.pdf")


def smooth_2d_distribution(x, y, window_size, num_points=100):
    """
    Smooth a 2D distribution of points using a co-moving window along the X-axis.

    Parameters:
    x (array-like): X-coordinates of the points.
    y (array-like): Y-coordinates of the points.
    window_size (float): Size of the moving window in X-axis units.
    num_points (int): Number of points in the smoothed output.

    Returns:
    tuple: (x_smooth, y_smooth) - Arrays of smoothed X and Y coordinates.
    """
    # Sort the points by x-coordinate
    sorted_indices = np.argsort(x)
    x_sorted = np.array(x)[sorted_indices]
    y_sorted = np.array(y)[sorted_indices]

    # Create evenly spaced points for the smoothed output
    x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), num_points)
    y_smooth = np.zeros_like(x_smooth)

    # Perform smoothing
    for i, x_i in enumerate(x_smooth):
        # Define window boundaries
        x_min = x_i - window_size / 2
        x_max = x_i + window_size / 2

        # Find points within the window
        mask = (x_sorted >= x_min) & (x_sorted <= x_max)

        # Calculate the mean y-value for points in the window
        if np.sum(mask) > 0:
            y_smooth[i] = np.mean(y_sorted[mask])
        else:
            # If no points in window, interpolate
            y_smooth[i] = np.interp(x_i, x_sorted, y_sorted)

    return x_smooth, y_smooth


def smooth_3d_distribution(x, y, z, window_size, num_points=100):
    # Sort the points by x-coordinate
    sorted_indices = np.argsort(x)
    x_sorted = np.array(x)[sorted_indices]
    y_sorted = np.array(y)[sorted_indices]
    z_sorted = np.array(z)[sorted_indices]

    # Create evenly spaced points for the smoothed output
    x_smooth = np.linspace(x_sorted.min(), x_sorted.max(), num_points)
    y_smooth = np.zeros_like(x_smooth)
    z_smooth = np.zeros_like(x_smooth)

    # Perform smoothing
    for i, x_i in enumerate(x_smooth):
        # Define window boundaries
        x_min = x_i - window_size / 2
        x_max = x_i + window_size / 2

        # Find points within the window
        mask = (x_sorted >= x_min) & (x_sorted <= x_max)

        # Calculate the mean y-value for points in the window
        print(y_sorted[mask])
        if np.sum(mask) > 0:
            y_smooth[i] = np.mean(y_sorted[mask])
            z_smooth[i] = np.mean(z_sorted[mask])
        else:
            # If no points in window, interpolate
            y_smooth[i] = np.interp(x_i, x_sorted, y_sorted)
            z_smooth[i] = np.interp(x_i, x_sorted, z_sorted)

    return x_smooth, y_smooth, z_smooth


def binned_mean_data(r, x):
    R = np.arange(0, r[-1], 0.1)
    X = np.zeros(len(R))
    N = np.zeros(len(R))
    for i in range(len(R) - 1):
        for j in range(len(r)):
            if r[j] >= R[i] and r[j] <= R[i + 1]:
                X[i] += x[j]
                N[i] += 1.0
    for i in range(len(X)):
        if X[i] > 0 and N[i] > 0:
            X[i] = X[i] / float(N[i])

    return R, X


def plot_ionization_fraction(pos, xion, temperature):
    T = [] | units.K
    r = [] | units.parsec
    x = []
    for pi, xi, Ti in zip(pos, xion, temperature):
        T.append(Ti)
        r.append(pi.length())
        x.append(xi)
    r, x, T = list(zip(*sorted(zip(r.value_in(units.parsec), x, T.value_in(units.K)))))

    R, X, Z = smooth_3d_distribution(r, x, T, window_size=0.04, num_points=20)
    print(R, X)

    figure = plt.figure(figsize=(10, 6))
    plt.rcParams.update({"font.size": 18})
    ax = figure.add_subplot(111)
    fig_fake = plt.figure()
    ax_fake = fig_fake.add_subplot(111)
    ax.set_xlabel("r [pc]")
    ax.set_ylabel(r"$\xi_{\rm ion}$")
    cbar = ax.scatter(r, x, c=T, lw=0, s=100, cmap="rainbow", alpha=0.1)

    figure.colorbar(cbar, extend="neither", label="T [K]")

    ax.plot(R, X, lw=2)
    ax.set_xlim((0, 7))
    ax.set_ylim((-0.04, 1.19))
    figure.savefig("fig_ionization_of_GMC.pdf")


# #BOOKLISTSTART1# #
def generate_ism_initial_conditions(number_of_particles, boxsize):
    converter = nbody_system.nbody_to_si(1.0 | units.MSun, 0.5 | units.parsec)
    ism = new_plummer_gas_model(number_of_particles, converter)

    mass_disk = 0.2 * ism.mass.sum()
    n_disk = int(mass_disk / ism[0].mass)
    converter = nbody_system.nbody_to_si(mass_disk, 1 | units.parsec)

    from amuse.ext.protodisk import ProtoPlanetaryDisk

    disk = ProtoPlanetaryDisk(
        n_disk,
        convert_nbody=converter,
        densitypower=1.5,
        Rmin=1,
        Rmax=5,
        q_out=2.0,
        discfraction=0.1,
    ).result
    disk.x += 1 | units.pc
    ism.add_particles(disk)
    ism.xion = 0.5

    hydro = Fi(converter)
    hydro.gas_particles.add_particles(ism)
    hydro.evolve_model(1 | units.hour)
    hydro.gas_particles.new_channel_to(ism).copy()
    hydro.stop()
    ism = ism.select(lambda r: r.length() < 0.5 * boxsize, ["position"])
    print(
        f"Max density: {ism.rho.max().in_(units.MSun/units.parsec**3)} ",
        f"Mean density: {np.mean(ism.rho.value_in(units.g/units.cm**3))} g/cc ",
        f"{ism.rho.max().in_(units.amu/units.cm**3)}",
        f"{np.mean(ism.rho.value_in(units.amu/units.cm**3))} amu/cc",
    )
    return ism, converter


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART# #
def rad_minimal(
    number_of_particles=10000,
    luminosity=100 | units.LSun,
    boxsize=100 | units.parsec,
    time_end=0.1 | units.Myr,
    seed=12345,
):
    np.random.seed(seed)
    ism, converter = generate_ism_initial_conditions(number_of_particles, boxsize)

    source = Particle()
    source.position = (1, 0, 0) | units.parsec
    source.luminosity = luminosity / (20.0 | units.eV)

    radiative = Sphray(convert=converter, mode="openmp")
    radiative.set_boundary(0)  # vacuum
    radiative.set_spectra_file("spectra/thermal6e4.cdf")
    radiative.set_raynumber(1.0e4 | units.Myr**-1)
    print(radiative.parameters)

    radiative.src_particles.add_particle(source)
    radiative.gas_particles.add_particles(ism)

    radiative.evolve_model(time_end)

    print(radiative.gas_particles)
    mux = mu()
    temperature = mux / constants.kB * radiative.gas_particles.u

    print(f"min ionization: {radiative.gas_particles.xion.min()}")
    print(f"average ionization: {radiative.gas_particles.xion.mean()}")
    print(f"max ionization: {radiative.gas_particles.xion.max()}")
    plot_ionization_fraction(
        radiative.gas_particles.position, radiative.gas_particles.xion, temperature
    )

    plot_gas_temperature(radiative)

    radiative.stop()


# #BOOKLISTSTOP# #


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-N",
        "--number_of_particles",
        type=int,
        default=10000,
        help="number of gas particles",
    )
    parser.add_argument(
        "-t",
        "--time_end",
        type=units.Myr,
        default=0.1 | units.Myr,
        help="radiation time",
    )
    parser.add_argument(
        "-L",
        "--luminosity",
        type=units.LSun,
        default=100 | units.LSun,
        help="luminosity of ionizing source",
    )
    parser.add_argument(
        "-d",
        "--boxsize",
        type=units.parsec,
        default=100 | units.parsec,
        help="size of the density box",
    )
    parser.add_argument("-s", "--seed", type=int, default=12345, help="random seed")

    return parser


def main(**kwargs):
    rad_minimal(**kwargs)


if __name__ == "__main__":
    arguments = new_argument_parser().parse_args()
    main(**arguments.__dict__)
