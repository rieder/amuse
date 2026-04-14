"""
Example to plot the performance of different Nbody integrators
"""
import argparse
import math
import matplotlib.pyplot as plt

lstyles = ['-.', '-.', '-.', '-.', '--', '-', '-', '-', '-']
lwidth = [2, 2, 4, 2, 4, 2, 2, 2, 4]

def read_file(filename, column, keyword):
    x = []
    with open(filename, encoding="utf-8") as fptr:
        lines = fptr.readlines()
    for line in lines:
        line = line.split()
        if line[1] == keyword:
            # print line
            x.append(float(line[column]))
    return x

def plot_nbody_performance(filename=None, lim=-1):
    integrators = [
        "Hermite", "MI6", "ph4", "Huayno", "ph4_GPU", "Gadget2", "BHTree", "Bonsai"
    ]

    if filename is None:
        return

    x_label = "N"
    y_label = "$t_{wall} [s]$"
    figure = plt.figure(figsize=(14, 10))
    ax = figure.add_subplot(1, 1, 1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    npp = [12, 512]
    ntc = [1024, 2000*1024]
    tpp = [0.01*x*x for x in npp]
    ttc = [1.e-6*(x*math.log(x)) for x in ntc]
    ax.plot(npp, tpp, c='k', ls='-', lw=4)
    ax.plot(ntc, ttc, c='k', ls='-', lw=4)
    ax.text(12, 8, '$N^2$')
    ax.text(2.e+3, 0.005, '$N \log (N)$')

    for ii, integrator in enumerate(integrators):

        x = read_file(filename, 6, integrator)
        y1 = read_file(filename, 22, integrator)
        y2 = read_file(filename, 23, integrator)
        for i in range(len(y1)):
            y1[i] *= 0.1
            y2[i] *= 0.1
        if len(x) > 0:
            ax.plot(x, y2, label=integrator, lw=lwidth[ii], ls=lstyles[ii])
            ax.scatter(x, y2, lw=0, s=200)

    ax.legend(loc="lower right", fontsize=18)

    save_file = "Nbody_performance.png"
    plt.savefig(save_file)
    print("\nSaved figure in file", save_file, '\n')
    plt.show()


def new_argument_parser():
    result = argparse.ArgumentParser()
    result.add_argument(
        "-f", dest="filename", default='NbodyAMUSE.test',
        help="output filename [NbodyAMUSE.test]"
        )
    result.add_argument(
        "-l", dest="lim", type=float, default=-1, help="boxsize"
    )
    return result


def main():
    arguments = new_argument_parser().parse_args()
    plot_nbody_performance(**arguments.__dict__)


if __name__ == "__main__":
    main()
