import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units
from amuse.support.console import set_printing_strategy
from amuse.datamodel import Particles, Particle
from amuse.community.seba import Seba

set_printing_strategy(
    "custom",  # nbody_converter = converter,
    preferred_units=[units.MSun, units.au, units.Myr],
    precision=4,
    prefix="",
    separator=" [",
    suffix="]",
)


def maximum_stellar_radius():
    M = (10 ** np.arange(-0.2, 2, 0.05)) | units.MSun
    R = [] | units.au
    for mi in M:
        stellar = Seba()
        s = Particle()
        s.mass = mi
        stellar.particles.add_particle(s)
        rmax = 0 | units.au
        while stellar.particles[0].stellar_type < 9 | units.stellar_type:
            stellar.evolve_model()
            rmax = max(rmax, stellar.particles[0].radius)
        print(
            "Time=",
            mi,
            stellar.particles[0].mass,
            stellar.model_time,
            stellar.particles[0].radius,
        )
        R.append(rmax)
    return M, R


names = [
    r"$\tau$ CMa",
    "V* CQ Dra",
    r"$\xi$ Tau",
    "V* V1334 Cyg",
    "V* DL Vir",
    "V* d Ser",
    "HD 97131",
    "KIC002856960",
]
Mprim = [50.0, 6.3, 5.5, 4.4, 3.68, 2.5, 1.5, 0.76] | units.MSun
Rprim = [0.95, 2.64, 0.42, 3.99, 2.34, 0.62, 0.26, 0.12] | units.au

figure = plt.figure()
ax = figure.add_subplot(1, 1, 1)
ax.minorticks_on()  # switch on the minor ticks
ax.locator_params(nbins=3)

M, Rmax = maximum_stellar_radius()

x_label = "m [$M_\odot$]"
y_label = "$r [au]$"
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.semilogx()
plt.semilogy()
plt.xlim(0.6, 100.0)
plt.ylim(0.07, 20.0)
plt.scatter(
    Mprim.value_in(units.MSun), Rprim.value_in(units.au), lw=0, s=300
)
plt.plot(M.value_in(units.MSun), Rmax.value_in(units.au), lw=4)
for i in range(len(names)):
    plt.text(
        1.1 * Mprim[i].value_in(units.MSun),
        0.9 * Rprim[i].value_in(units.au),
        names[i],
        fontsize=20,
    )

save_file = "fig_RLOF_triples_R.png"
plt.savefig(save_file)
print("\nSaved figure in file", save_file, "\n")
plt.show()
