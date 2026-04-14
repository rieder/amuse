"""
Gravity drifter example
"""
import numpy as np
import matplotlib.pyplot as plt
from amuse.datamodel import Particles
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.ic.kingmodel import new_king_model


# #BOOKLISTSTART1# #
class drift_without_gravity:
    def __init__(self, convert_nbody, time=0 | units.Myr):
        self.model_time = time
        self.convert_nbody = convert_nbody
        self.particles = Particles()

    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity * dt
        self.model_time = t_end

    @property
    def potential_energy(self):
        return quantities.zero

    @property
    def kinetic_energy(self):
        return (
            0.5 * self.particles.mass * self.particles.velocity.lengths() ** 2
        ).sum()

    def stop(self):
        pass
# #BOOKLISTSTOP1# #


class MilkyWayGalaxy:

    def __init__(
        self,
        Mb=1.40592e10 | units.MSun,
        Md=8.5608e10 | units.MSun,
        Mh=1.07068e11 | units.MSun,
    ):
        self.Mb = Mb
        self.Md = Md
        self.Mh = Mh

    def get_potential_at_point(self, eps, x, y, z):
        """
        Returns gravitational potential at given location.
        `eps` is ignored, but kept as an argument for compatibility reasons.
        """
        r = (x**2 + y**2 + z**2) ** 0.5
        R = (x**2 + y**2) ** 0.5
        # bulge
        b1 = 0.3873 | units.kpc
        pot_bulge = -constants.G * self.Mb / (r**2 + b1**2) ** 0.5
        # disk
        a2 = 5.31 | units.kpc
        b2 = 0.25 | units.kpc
        pot_disk = (
            -constants.G * self.Md
            / (R**2 + (a2 + (z**2 + b2**2) ** 0.5) ** 2) ** 0.5
        )
        # halo
        a3 = 12.0 | units.kpc
        cut_off = 100 | units.kpc
        d1 = r / a3
        c = 1 + (cut_off / a3) ** 1.02
        pot_halo = -constants.G * (self.Mh / a3) * d1**1.02 / (1 + d1**1.02) - (
            constants.G * self.Mh / (1.02 * a3)
        ) * (
            -1.02 / c + np.log(c) + 1.02 / (1 + d1**1.02)
            - np.log(1.0 + d1**1.02)
        )
        return 2 * (pot_bulge + pot_disk + pot_halo)  # multiply by 2 for
        # a rigid potential

    def get_gravity_at_point(self, eps, x, y, z):
        """
        Returns gravitational acceleration at given location.
        `eps` is ignored, but kept as an argument for compatibility reasons.
        """
        r = (x**2 + y**2 + z**2) ** 0.5
        R = (x**2 + y**2) ** 0.5
        # bulge
        b1 = 0.3873 | units.kpc
        force_bulge = -constants.G * self.Mb / (r**2 + b1**2) ** 1.5
        # disk
        a2 = 5.31 | units.kpc
        b2 = 0.25 | units.kpc
        d = a2 + (z**2 + b2**2) ** 0.5
        force_disk = -constants.G * self.Md / (R**2 + d**2) ** 1.5
        # halo
        a3 = 12.0 | units.kpc
        d1 = r / a3
        force_halo = -constants.G * self.Mh * d1**0.02 / (a3**2 * (1 + d1**1.02))

        ax = force_bulge * x + force_disk * x + force_halo * x / r
        ay = force_bulge * y + force_disk * y + force_halo * y / r
        az = (
            force_bulge * z
            + force_disk * d * z / (z**2 + b2**2) ** 0.5
            + force_halo * z / r
        )

        return ax, ay, az

    def circular_velocity(self, r):
        z = 0 | units.kpc
        b1 = 0.3873 | units.kpc
        a2 = 5.31 | units.kpc
        b2 = 0.25 | units.kpc
        a3 = 12.0 | units.kpc

        rdphi_b = constants.G * self.Mb * r**2 / (r**2 + b1**2) ** 1.5
        rdphi_d = (
            constants.G
            * self.Md
            * r**2
            / (r**2 + (a2 + (z**2 + b2**2) ** 0.5) ** 2) ** 1.5
        )
        rdphi_h = (
            constants.G
            * self.Mh
            * (r / a3) ** 0.02
            * r
            / (a3**2 * (1 + (r / a3) ** 1.02))
        )

        vel_circb = rdphi_b
        vel_circd = rdphi_d
        vel_circh = rdphi_h

        return (vel_circb + vel_circd + vel_circh) ** 0.5


def plot_cluster(x, y, width=20 | units.kpc):
    """
    Simple function to plot a star cluster.
    """
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_xlabel("X [kpc]")
    ax.set_ylabel("Y [kpc]")
    ax.set_xlim(-0.5*width.value_in(units.kpc), 0.5*width.value_in(units.kpc))
    ax.set_ylim(-0.5*width.value_in(units.kpc), 0.5*width.value_in(units.kpc))

    plt.scatter(x, y, s=50, lw=0)
    plt.show()


def evolve_cluster_in_galaxy(
    number_of_stars,
    W0=3,
    radius_init=8.5 | units.kpc,
    time_end=4.6 | units.Gyr,
    time_step=0.1 | units.Myr,
    mass_cluster=5.0e4 | units.MSun,
    radius_cluster=1.0 | units.parsec,
):

    galaxy_code = MilkyWayGalaxy()

    converter = nbody_system.nbody_to_si(mass_cluster, radius_cluster)
    cluster_code = drift_without_gravity(convert_nbody=converter)
    bodies = new_king_model(number_of_stars, W0, convert_nbody=converter)
    cluster_code.particles.add_particles(bodies)

    stars = cluster_code.particles.copy()
    stars.x += radius_init
    stars.vy = 0.8 * galaxy_code.circular_velocity(radius_init)
    channel = stars.new_channel_to(cluster_code.particles)
    channel.copy_attributes(["x", "y", "z", "vx", "vy", "vz"])

    system = bridge(verbose=False)
    system.add_system(cluster_code, (galaxy_code,))

    times = quantities.arange(0 | units.Myr, time_end, 100 * time_step)
    for i, time in enumerate(times):
        print(f"Time={time.in_(units.Myr)}")
        system.evolve_model(time, timestep=time_step)

    x = system.particles.x.value_in(units.kpc)
    y = system.particles.y.value_in(units.kpc)
    cluster_code.stop()
    return x, y


def main():
    number_of_stars = 1024
    W0 = 3
    radius_init = 8.5 | units.kpc
    time_step = 0.1 | units.Myr
    time_end = 4.6 | units.Gyr
    mass_cluster = 5.0e4 | units.MSun
    radius_cluster = 1.0 | units.parsec
    x, y = evolve_cluster_in_galaxy(
        number_of_stars, W0=W0, radius_init=radius_init, time_end=time_end,
        time_step=time_step, mass_cluster=mass_cluster,
        radius_cluster=radius_cluster
    )
    plot_cluster(x, y)


if __name__ == "__main__":
    main()
