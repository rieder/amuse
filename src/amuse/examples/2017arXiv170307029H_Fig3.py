"""
Example code to reproduce Figure 3 of the paper 2017arXiv170307029H
"""

import argparse
import os
import matplotlib.pyplot as plt
from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.io import read_set_from_file


def plot_radial_distribution(gas, nbin, c, lw):
    X = []
    Ymean = []
    # Ystd = []
    for gi in range(len(gas)-nbin):
        X.append(gas[gi: gi+nbin].r.value_in(units.parsec).mean())
        S = gas[gi+nbin].r**2-gas[gi].r**2
        rho = gas[gi: gi+nbin].mass.sum()/S
        # if hasattr(gas, "rho"):
        # rho = gas[gi: gi+nbin].rho.max()*S.sqrt()
        Ymean.append(rho.value_in(units.MSun/units.parsec**2))
    plt.plot(X, Ymean, lw=lw)


def plot_radial_density_distribution(gas, stars):
    com = stars.center_of_mass()

    gas.r = ((gas.x-com[0])**2 + (gas.y-com[1])**2).sqrt()
    gas = gas.sorted_by_attributes("r")

    stars.r = ((stars.x-com[0])**2 + (stars.y-com[1])**2).sqrt()
    stars = stars.sorted_by_attributes("r")

    max_age = stars.birth_age.max()
    young_stars = stars[max_age-stars.birth_age < 1.5 | units.Myr]

    # mN2H = 2*29.02134
    # cutoff_density = (
    #     1e5 * mN2H * constants.atomic_unit_of_mass / (1 | units.cm**3)
    # )
    # dense_gas = gas[gas.density > cutoff_density]

    plot_radial_distribution(young_stars, 10, c=6, lw=4)
    plot_radial_distribution(stars, 60, c=0, lw=4)
    plot_radial_distribution(gas, 100, c=3, lw=4)
    # plot_radial_distribution(dense_gas, 100, c=7, lw=2)
    plt.xlim(0, 1.1)
    plt.ylim(1, 1100)
    plt.semilogy()
    plt.xlabel("R [pc]")
    plt.ylabel(r"$\Sigma$ [M$_\odot$ pc$^{-2}]$")
    plt.savefig("2017arXiv170307029H_Fig3")


def reproduce_2017arXiv170307029H_Fig3(filename=None):
    R_young_stars = [0.1012221753955787, 0.2968660204589863,
                     0.5017430508551541, 0.7046418772801384,
                     0.8984687805821738, 1.1087212228990806, 1.296848150063784]
    S_young_stars = [69.68021137071266, 21.28302762211248,
                     8.871846723535821, 4.634361707743547, 0.842276996426784,
                     0.35283784451473543, 0.32284672967371486]

    R_all_stars = [0.09966839700574376, 0.29542608638995443,
                   0.4910636790322119, 0.6989696157755654, 0.9093883458890231,
                   1.1027718150669223, 1.3041555749011087]

    S_all_stars = [246.4266804633498, 108.52660553059054,
                   59.62962877617595, 22.170446486030563,
                   6.83698178212555, 2.8123977037860786,
                   0.28189901519958876]

    R_dense_gas = [0.10160842129883879, 0.30957902551862887,
                   0.5005420834356191, 0.6940641257773106,
                   0.8993648870413651, 1.107381682315753,
                   1.2954901330586175]

    S_dense_gas = [587.057496278814, 193.757708949030, 81.05652212929158,
                   25.83103135675521, 13.72720495792530,
                   4.161094960653975, 3.939217817754274]

    R_all_gas = [0.09890162549942522, 0.2971536318319011,
                 0.5047269929821525, 0.6955052866807527,
                 0.8934247174201266, 1.0960741121502846,
                 1.296349286674131]

    S_all_gas = [1011.82389615108, 382.3950057709804, 262.3523868909657,
                 154.2500699106347,
                 107.5704842567941,
                 88.9640154049995,
                 68.72257895441606]

    plt.rcParams.update({'font.size': 30})
    fig, ax = plt.subplots(figsize=[16, 10])
    ax.minorticks_on()  # switch on the minor ticks
    ax.tick_params('both', length=15, width=2, which='major')
    ax.tick_params('both', length=6, width=1, which='minor')
    # colors = get_distinct(10)

    plt.scatter(R_young_stars, S_young_stars, s=100, marker='s', lw=0)
    plt.scatter(R_all_stars, S_all_stars, s=100, marker='s', lw=0)
    plt.scatter(R_all_gas, S_all_gas, s=100, marker='s', lw=0)
    # plt.scatter(R_dense_gas, S_dense_gas, s=100, marker='s', lw=0)
    disk = Particles(0)
    stars = Particles(0)

    if not os.path.exists(filename):
        # try to get the file from the data directory in the examples
        filename = f"data/{filename}"
        filename = os.path.join(os.path.dirname(__file__), filename)
    bodies = read_set_from_file(filename, "amuse")
    for bi in bodies.history:
        # print(bi)
        if len(bi) > 0:
            if hasattr(bi, "name") and "gas" in str(bi.name):
                disk.add_particles(bi.copy())
            elif "Star" in str(bi.name):
                stars.add_particles(bi.copy())

    print(
        "Stellar masses:",
        stars.mass.min().in_(units.MSun),
        stars.mass.mean().in_(units.MSun),
        stars.mass.max().in_(units.MSun),
        stars.mass.median().in_(units.MSun),
    )

    # plot_age_gasdensity(disk, stars)
    plot_radial_density_distribution(disk, stars)


def new_argument_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-f",
        "--filename",
        default="GMC_R2pcN20k_SE_T45Myr.amuse",
        help="output filename",
    )
    return parser


def main():
    args = new_argument_parser().parse_args()
    reproduce_2017arXiv170307029H_Fig3(**args.__dict__)


if __name__ == "__main__":
    main()
