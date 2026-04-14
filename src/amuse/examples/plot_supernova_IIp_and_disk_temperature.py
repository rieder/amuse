import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units

from amuse.examples.supernova_IIp_Lightcurve import Supernova_IIp


def read_supernova_irradiation_file(filename):
    t = []
    Tmean = []
    Tmax = []
    for line in open(filename):
        if "Time" in line:
            sl = line.split()
            t.append(float(sl[1]))
        if "Temperature:" in line:
            sl = line.split()
            Tmean.append(float(sl[3]))
            Tmax.append(float(sl[1]))
    return np.asarray(t, dtype="float"), Tmean, Tmax


def read_Earthorbit():
    t = []
    e = []
    a = []
    tmax = 1.0e6
    for line in open("Eart_Orbit_Eps-3.data"):
        if "#" not in line:
            sl = line.split()
            t.append(float(sl[0]))
            a.append(float(sl[1]))
            e.append(float(sl[2]))
            if t[-1] > tmax:
                break
    return t, a, e


def plot_supernova_IIp_and_disk_temperature():
    to = 50 | units.day

    t_offset = to + (((0.15 | units.parsec) / (1 | units.lightyear)) | units.yr)
    filename = "SN10a.R0.15.i15.data"
    time, Tmean, Tmax = read_supernova_irradiation_file(filename)
    time += t_offset.value_in(units.day)

    t_offset = to + (((0.3 | units.parsec) / (1 | units.lightyear)) | units.yr)
    filename = "SN11aof.R0.3.i45.data"
    t3pc_N7, Tmean3pc_N7, Tmax3pc_N7 = read_supernova_irradiation_file(filename)
    t3pc_N7 += t_offset.value_in(units.day)

    t_offset = to + (((0.4 | units.parsec) / (1 | units.lightyear)) | units.yr)
    filename = "SN11aof.R0.4.i15.data"
    t3pc_N8, Tmean3pc_N8, Tmax3pc_N8 = read_supernova_irradiation_file(filename)
    t3pc_N8 += t_offset.value_in(units.day)

    PS1_11aof = Supernova_IIp("11aof", to)
    t = 10 ** np.arange(-2, 3.0, 0.01) | units.day
    L11aof = [] | units.erg / units.s
    for ti in t:
        L11aof.append(PS1_11aof.luminosity_at_time(ti))
    L11aof = np.log10(L11aof.value_in(units.LSun))

    PS1_10a = Supernova_IIp("10a", to)
    L10a = [] | units.erg / units.s
    for ti in t:
        L10a.append(PS1_10a.luminosity_at_time(ti))
    L10a = np.log10(L10a.value_in(units.LSun))

    x_label = "$t$ [day]"
    y_label = "L [L$_\odot$]"
    figure = plt.figure()
    ax1 = figure.add_subplot(1, 1, 1)
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ax1.plot(t.value_in(units.day), L11aof, ls="-", c=colors[0])
    ax1.plot(t.value_in(units.day), L10a, ls="--", c=colors[0])
    ax1.set_xlabel("time [day]")
    ax1.set_ylabel("log$_{10}$(L/L$_\odot$)", color=colors[0])
    for tl in ax1.get_yticklabels():
        tl.set_color(colors[0])

    ax2 = ax1.twinx()
    ax2.plot(time, Tmean, colors[1], ls="--")
    ax2.plot(t3pc_N7, Tmean3pc_N7, colors[1])
    ax2.plot(t3pc_N8, Tmean3pc_N8, colors[1], lw=4)
    ax2.set_ylabel("mean temperature [K]", color=colors[1])
    for tl in ax2.get_yticklabels():
        tl.set_color(colors[1])

    t_cooling = [950, 1061]
    T_cooling = [1600, 800]
    ax2.plot(t_cooling, T_cooling, colors[3], lw=1)
    ax2.text(
        t_cooling[0] + 20,
        T_cooling[0] - 100,
        "cooling of 0.3 K/h",
        rotation=-76.5,
        color=colors[3],
    )

    plt.savefig("supernova_IIp_and_disk_temperature")


def main():
    plot_supernova_IIp_and_disk_temperature()


if __name__ == "__main__":
    main()
