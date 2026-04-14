from amuse.lab import *
bodies = read_set_from_file("../data/Pl_N16NO_A1000R.amuse", close_file=True)
stars = bodies[bodies.name=="star"]
systems = bodies[bodies.name=="system"]
planets = bodies[bodies.name=="planet"]
debris = bodies[bodies.name=="asteroid"]

print(f"N={len(bodies)}, {len(stars)}, {len(systems)}, {len(planets)}, {len(debris)}")

from matplotlib import pyplot as plt

figure = plt.figure(figsize=(6, 6))
ax = plt.gca()
ax.minorticks_on() # switch on the minor ticks
ax.locator_params(nbins=3)

plt.xlabel("x [pc]")
plt.ylabel("y [pc]")
plt.xlim(-0.01, 0.01)
plt.ylim(-0.01, 0.01)
plt.scatter(debris.x.value_in(units.pc), debris.y.value_in(units.pc), c='k', lw=0, s=1, label="debris")
plt.scatter(planets.x.value_in(units.pc), planets.y.value_in(units.pc), c='y', lw=0, s=20, label="planet")
plt.scatter(systems.x.value_in(units.pc), systems.y.value_in(units.pc), c='r', lw=0, s=30, label="planetary system")
plt.scatter(stars.x.value_in(units.pc), stars.y.value_in(units.pc), c='b', lw=0, s=30, label="single star")
plt.savefig("fig_Pl_N16NO_A1000R.pdf")
plt.show()
