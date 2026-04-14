import os
import pickle
import matplotlib.pyplot as plt


def main():
    datadir = os.path.dirname(__file__)
    with open(
        os.path.join(datadir, "data") + "Obs_Trapezium_disksizes.pkl",
        "r",
        encoding="utf-8",
    ) as filepointer:
        (radius_obs, yc_obs) = pickle.load(filepointer)

    with open(
        os.path.join(datadir, "data") + "Tr_N2000_R0.5pc_Q0.5_F1.6.pkl",
        "r",
        encoding="utf-8",
    ) as filepointer:
        (radius_sim, yc_sim) = pickle.load(filepointer)

    print(len(radius_obs), len(yc_obs))
    print(len(radius_sim), len(yc_sim))

    color = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    x_label = r"R [$R_\odot$]"
    y_label = "$f_{<R}$"
    figure = plt.figure(figsize=(8, 8))
    ax = figure.add_subplot(111)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    ax.plot(radius_obs, 95 * yc_obs, c=color[0])
    ax.plot(radius_sim, 95 * yc_sim, c=color[1], ls="--")

    plt.savefig("Tr_N2000_R05pc_Q05_F16_r1")


if __name__ == "__main__":
    main()
