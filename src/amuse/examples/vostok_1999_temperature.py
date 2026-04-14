import sys
import numpy as np
import matplotlib.pyplot as plt
if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files


def read_vostok_1999_temperature_data(datafile):
    t = []
    T = []
    for line in open(datafile):
        if "***" not in line and "#" not in line:
            sl = line.split()
            t.append(float(sl[1]))
            T.append(float(sl[3]))
    return t, T


def read_earth_orbit(datafile):
    t = []
    e = []
    a = []
    tmax = 1.0e6
    for line in open(datafile):
        if "#" not in line:
            sl = line.split()
            t.append(float(sl[0]))
            a.append(float(sl[1]))
            e.append(float(sl[2]))
            if t[-1] > tmax:
                break
    return t, a, e


def vostok_1999_temperature():
    data_dir = files("amuse_examples.data")
    tv, T = read_vostok_1999_temperature_data(
        f"{data_dir}/vostok_1999_temperature.data"
    )
    te, ae, ee = read_earth_orbit(f"{data_dir}/earth_orbit_eps_3.data")

    q = []
    q0 = ae[0] * (1 - ee[0] ** 2)
    ao = ae[0]
    T_mean = np.mean(T)
    Tmin = np.min(T)
    Tmax = np.max(T)
    dT = Tmax - Tmin
    for i in range(len(te)):
        te[i] -= 0
        ae[i] = 2e5 * (ae[i] - ao) / ao
        te[i] = te[i] / 1000.0
    for i in range(len(tv)):
        tv[i] = tv[i] / 1000.0

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    figure = plt.figure()
    ax1 = figure.add_subplot(1, 1, 1)
    ax1.set_xlabel("$a-a_0$ [au]")
    ax1.set_ylabel("eccentricty")

    ax1.plot(tv, T, ls="-", c=colors[0])
    ax1.set_xlabel("time before present (kyr)")
    ax1.set_xlim((450, 0))

    # Make the y-axis label and tick labels match the line color.

    ax1.set_ylabel("temperature variation [$^o$C]", color=colors[0])
    for tl in ax1.get_yticklabels():
        tl.set_color(colors[0])

    ax2 = ax1.twinx()
    ax2.plot(te, ee, colors[1])
    ax2.set_ylabel("eccentricity", color=colors[1])
    for tl in ax2.get_yticklabels():
        tl.set_color(colors[1])

    save_file = "vostok_1999_temperature.pdf"
    plt.savefig(save_file)
    print(f"\nSaved figure in file {save_file}\n")


def main(**kwargs):
    vostok_1999_temperature(**kwargs)


if __name__ == "__main__":
    main()
