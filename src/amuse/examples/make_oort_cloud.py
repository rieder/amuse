"""
Spherical isotropic (Oort) cloud -- initial conditions for a particle set

See for example:
  * Duncan, M.; Quinn, T.; Tremaine, S.
  -- http://adsabs.harvard.edu/abs/1987AJ.....94.1330D
  * Feng, F.; Bailer-Jones, C. A. L. 
  -- http://adsabs.harvard.edu/abs/2014MNRAS.442.3653F
  * Dybczynski, P. A. -- http://adsabs.harvard.edu/abs/2002A%26A...396..283D

Initial conditions -- given by distributions of orbital elements:
  semi-major axes -- power law, default: dn/da ~ a^(-1.5),
  eccentricities -- dn/de ~ e,
  constrain on the minimum pericenter,
  isotropic orbits -- distribution of orbital inclinations: cos(i) = -1--1,
      longitude of ascending node: 0--2pi,
      argument of periastron: 0--2pi,
      mean anomaly: 0--2pi,
  equal mass particles.
"""
import argparse
import numpy as np

from amuse.units import units, nbody_system
from amuse.datamodel import Particles
from amuse.community.kepler import Kepler
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.io import read_set_from_file, write_set_to_file


def random_power_min_max(size, x_min, x_max, exp_plus_one):
    """
    returns random floats in the interval [x_min,x_max] drawn from distribution
    pdf(x) = const * x**(exp_plus_one-1), x_min <= x <= x_max;
    assuming: x_min < x_max, exp_plus_one != 0
    """
    r = np.random.random(size=size)
    x_min_gamma = x_min**exp_plus_one
    x_max_gamma = x_max**exp_plus_one
    return (x_min_gamma + (x_max_gamma - x_min_gamma) * r) ** (1.0 / exp_plus_one)


def relative_position_and_velocity_from_orbital_elements(
    mass1, mass2, semimajor_axis, eccentricity, mean_anomaly, seed=None
):
    """
    Function that returns relative positions and velocity vectors or orbiters
    with masses mass2 of the central body with mass mass1 in Cartesian
    coordinates; for vectors of orbital elements -- semi-major axes,
    eccentricities, mean anomalies.  3D orientation of orbits (inclination,
    longitude of ascending node and argument of periapsis) are random.
    (cos(incl) is uniform -1--1, longitude of ascending node and argument of
    periapsis are uniform 0--2pi) Assuming mass1 is static in the center
    [0,0,0] m, [0,0,0] km/s (that is mass2<<mass1)
    """
    position_vectors = []
    velocity_vectors = []
    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.au)
    kepler = Kepler(converter)
    kepler.initialize_code()
    r_vec = (0.0, 0.0, 0.0) | units.au
    v_vec = (0.0, 0.0, 0.0) | units.kms
    # to change seed for each particle
    if seed is not None:
        i = 0
    for m2_i, a_i, ecc_i, ma_i in zip(
        mass2, semimajor_axis, eccentricity, mean_anomaly
    ):
        # print m2_i, a_i, ecc_i, ma_i
        if seed is not None:
            kepler.set_random(seed + i)
            i = i + 1
        kepler.initialize_from_elements(
            mass=(mass1 + m2_i),
            semi=a_i,
            ecc=ecc_i,
            mean_anomaly=ma_i,
            random_orientation=-1,
        )
        ri = kepler.get_separation_vector()
        vi = kepler.get_velocity_vector()
        # this is to get ~half of the orbits retrograde (that is with inclination
        # of 90--180 degrees) --> velocity = -velocity
        vel_vec_dir = np.random.random()
        if vel_vec_dir <= 0.5:
            vel_orientation = 1.0
        else:
            vel_orientation = -1.0
        position_vectors.append([ri[0], ri[1], ri[2]])
        velocity_vectors.append(
            [vel_orientation * vi[0], vel_orientation *
                vi[1], vel_orientation * vi[2]]
        )
    kepler.stop()
    return position_vectors, velocity_vectors


def ecc_random_power_with_min_peri(n, semi, min_peri, power=2.0):
    """
    random distribution in eccentricity P(e)~e
    (power = actual_exponent + 1)
    with minimum pericenter of min_peri for given semi-major axes semi
    """
    x = np.random.power(power, size=n)
    peri = semi * (1.0 - x)
    while np.any(peri < min_peri):
        filter_small_peri = peri < min_peri
        n_new = sum(filter_small_peri)
        # print "\t updating q", peri.min().in_(units.au), n_new
        x_random_new = np.random.power(power, size=n_new)
        x_new = 1.0 * x
        x_new[filter_small_peri] = x_random_new
        x = 1.0 * x_new
        peri = semi * (1.0 - x)
    return x


def get_cloud_orbital_elements(
    target_number_of_particles,
    m_total=0.0 | units.MSun,
    a_min=3000.0 | units.au,
    a_max=100000.0 | units.au,
    q_min=32.0 | units.au,
    gamma=-1.5,
):
    if q_min is not None:
        if q_min > a_min:
            a_min_q_corr = q_min
        else:
            a_min_q_corr = a_min
    else:
        a_min_q_corr = a_min

    a_in_au = random_power_min_max(
        target_number_of_particles, a_min_q_corr.value_in(
            units.au), a_max.value_in(units.au), gamma + 1.0
    )
    a = a_in_au * 1.0 | units.au
    ecc = ecc_random_power_with_min_peri(
        target_number_of_particles, a, q_min, power=2.0,
    )
    m_comets = (m_total / target_number_of_particles) * np.ones_like(ecc)
    return a, ecc, m_comets


class SphericalIsotropicCloudWithPlanets():
    def __init__(
        self,
        target_number_of_particles,
        m_star=1.0 | units.MSun,
        m_cloud=0.0 | units.MSun,
        a_min=3000.0 | units.au,
        a_max=100000.0 | units.au,
        q_min=32.0 | units.au,
        gamma=-1.5,
        seed=None,
    ):
        self.target_number_of_particles = target_number_of_particles
        self.m_star = m_star
        self.m_cloud = m_cloud
        self.a_min = a_min
        self.a_max = a_max
        self.q_min = q_min
        self.gamma = gamma
        self.seed = seed

        if self.q_min is not None:
            if self.q_min > self.a_min:
                self.a_min_q_corr = self.q_min
            else:
                self.a_min_q_corr = self.a_min
        else:
            self.a_min_q_corr = self.a_min

    def new_model(self):
        if self.seed is not None:
            np.random.seed(self.seed)

        a_in_au = random_power_min_max(
            self.target_number_of_particles,
            self.a_min_q_corr.value_in(units.au),
            self.a_max.value_in(units.au),
            self.gamma + 1.0,
        )
        a = a_in_au * 1.0 | units.au
        ecc = ecc_random_power_with_min_peri(
            self.target_number_of_particles, a, self.q_min, power=2.0)
        mean_anomaly = 2.0 * np.pi * np.random.random(
            size=self.target_number_of_particles,
        )
        m_comets = (self.m_cloud / self.target_number_of_particles) * np.ones_like(ecc)

        position_vectors, velocity_vectors = (
            relative_position_and_velocity_from_orbital_elements(
                self.m_star, m_comets, a, ecc, mean_anomaly, seed=self.seed
            )
        )
        return (m_comets, position_vectors, velocity_vectors)

    @property
    def result(self):
        masses, position_vectors, velocity_vectors = self.new_model()
        result = Particles(self.target_number_of_particles)
        result.mass = masses
        result.position = position_vectors
        result.velocity = velocity_vectors
        return result


def new_isotropic_cloud_with_planets(
    number_of_particles, *list_arguments, **keyword_arguments
):
    """
    Spherical isotropic cloud ~ Oort cloud given by distributions of orbital elements:
      semi-major axes -- power law, default: dn/da ~ a^(-1.5),
      eccentricities -- dn/de ~ e,
      constrain on the minimum pericenter,
      isotropic orbits -- distribution of orbital inclinations: cos(i) = -1--1,
          longitude of ascending node: 0--2pi,
          argument of periastron: 0--2pi,
          mean anomaly: 0--2pi,
      equal mass particles

    Cloud with planets : n_planets number of particles have m_planets masses
      -- next to equal mass comets, there is n_planets particles with masses m_planets
      -- total number of particles is number_of_particles
      -- total mass is: m_cloud - (m_cloud/number_of_particles)*n_planets + sum(m_planets)
          for number of_particles >> 1, total mass ~ m_cloud + sum(m_planets)

    The default values correspond to papers:
    * Duncan, M.; Quinn, T.; Tremaine, S. -- http://adsabs.harvard.edu/abs/1987AJ.....94.1330D
    * Feng, F.; Bailer-Jones, C. A. L. -- http://adsabs.harvard.edu/abs/2014MNRAS.442.3653F (their DQT model)

    :argument number_of_particles: number of particles to include in the cloud; including planets
    :argument m_star: mass of the central star
    :argument m_cloud: total mass of the cloud (particles are equal mass)
    :argument a_min: minimal semimajor axis
    :argument a_max: maximal semimajor axis, a_min < a_max
    :argument q_min: minimal pericenter
    :argument gamma: exponent of the semimajor axis distribution, f(a) ~ a^(gamma)
    :argument seed: random seed -- set to reproduce exactly the same IC
    """
    comets_and_planets = SphericalIsotropicCloudWithPlanets(
        number_of_particles, *list_arguments, **keyword_arguments
    ).result
    return comets_and_planets


def generate_initial_oort_cloud(m_star, n_comets, q_min, a_min, a_max, seed):
    """
    Generate a random Oort cloud
    save the result in Cartesian coordinates
    in default Rebound N-body units
    """
    cloud = SphericalIsotropicCloudWithPlanets(
        target_number_of_particles=n_comets,
        m_star=m_star,
        m_cloud=0.0 | units.MSun,
        a_min=a_min,
        a_max=a_max,
        q_min=q_min,
        gamma=-1.5,
        seed=seed,
    ).result

    r_time_unit = (1.0 | units.yr) / (2.0 * np.pi)

    comet_radius = 0.001 | units.RSun
    cloud_xyz = np.empty([n_comets, 8])
    cloud_xyz[:, 0] = cloud.mass.value_in(units.kg)
    cloud_xyz[:, 1] = comet_radius.value_in(units.au)
    cloud_xyz[:, 2] = cloud.x.value_in(units.au)
    cloud_xyz[:, 3] = cloud.y.value_in(units.au)
    cloud_xyz[:, 4] = cloud.z.value_in(units.au)
    cloud_xyz[:, 5] = (cloud.vx * r_time_unit).value_in(units.au)
    cloud_xyz[:, 6] = (cloud.vy * r_time_unit).value_in(units.au)
    cloud_xyz[:, 7] = (cloud.vz * r_time_unit).value_in(units.au)
    print(f" ** random cloud : M={m_star.in_(units.MSun)}"
          f"N {n_comets} seed {seed}")
    # print cloud_xyz[:10,:]

    return cloud_xyz


def hill_radius(mass_star, semi_major_axis, mass_planet):
    return semi_major_axis * (mass_planet / (3.0 * mass_star)) ** (1.0 / 3.0)


def add_comets(star, m_comets, n_comets, q_min, a_min, a_max, seed):
    # bodies.position -= star.position
    # bodies.velocity -= star.velocity

    cloud_xyz = generate_initial_oort_cloud(
        star.mass, n_comets, q_min, a_min, a_max, seed
    )
    masses = new_salpeter_mass_distribution(
        n_comets, 0.1 | units.MSun, 100 | units.MSun
    )
    masses = m_comets * masses / masses.sum()
    print("total mass in comets: M=", masses.sum().in_(units.MEarth))

    comets = Particles(n_comets)
    comets.mass = masses
    comets.name = "comet"
    comets.host = star
    comets.x = cloud_xyz[:, 2] | units.au
    comets.y = cloud_xyz[:, 3] | units.au
    comets.z = cloud_xyz[:, 4] | units.au
    comets.vx = cloud_xyz[:, 5] | 2 * np.pi * units.au / units.yr
    comets.vy = cloud_xyz[:, 6] | 2 * np.pi * units.au / units.yr
    comets.vz = cloud_xyz[:, 7] | 2 * np.pi * units.au / units.yr

    comets.position += star.position
    comets.velocity += star.velocity

    return comets


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument(
        "-f",
        dest="filename",
        default="starandplanets.amuse",
        help="input filename",
    )
    result.add_argument(
        "-F", dest="outfile", default=None, help="output filename"
    )
    result.add_argument(
        "--n_comets",
        dest="n_comets",
        type=int,
        default=100,
        help="number of comets",
    )
    result.add_argument(
        "--m_comets",
        type=units.MEarth,
        dest="m_comets",
        default=0 | units.MEarth,
        help="mass of all comets for a 1MSun star",
    )
    result.add_argument(
        "--seed",
        dest="seed",
        type=int,
        default=None,
        help="random number seed",
    )
    result.add_argument(
        "--q_min",
        type=units.au,
        dest="q_min",
        default=32 | units.au,
        help="initial minimal pericenter distance",
    )
    result.add_argument(
        "--a_min",
        type=units.au,
        dest="a_min",
        default=3000.0 | units.au,
        help="initial maximum distance",
    )
    result.add_argument(
        "--a_max",
        type=units.au,
        dest="a_max",
        default=100000.0 | units.au,
        help="initial maximum distance",
    )
    result.add_argument(
        "--m_min",
        type=units.MSun,
        dest="m_min",
        default=0.0 | units.MSun,
        help="maximum stellar mass with comets",
    )
    result.add_argument(
        "--m_max",
        type=units.MSun,
        dest="m_max",
        default=100.0 | units.MSun,
        help="maximum stellar mass with comets",
    )
    return result


def make_oort_cloud(
    outfile=None,
    filename="starandplanets.amuse",
    n_comets=100,
    m_comets=0 | units.MEarth,
    seed=None,
    q_min=32 | units.au,
    a_min=3000.0 | units.au,
    a_max=100000.0 | units.au,
    m_min=0.0 | units.MSun,
    m_max=100.0 | units.MSun,
):
    if outfile is None:
        outfile = f"sun_and_comets_N{n_comets:01}.amuse"

    bodies = read_set_from_file(filename, close_file=True)
    stars = bodies[bodies.type == "star"]
    selected_stars = stars.select(lambda m: m > m_min, ["mass"])
    selected_stars = selected_stars.select(lambda m: m <= m_max, ["mass"])
    print(f"Number of stars: {len(stars)} {len(selected_stars)}")

    for si in selected_stars:
        m_comets = m_comets * si.mass.value_in(units.MSun)
        n_comets = int(n_comets * si.mass.value_in(units.MSun))
        comets = add_comets(si, m_comets, n_comets,
                            q_min, a_min, a_max, seed)
        bodies.add_particles(comets)
    bodies.move_to_center()
    print(bodies)

    time = 0 | units.Myr
    write_set_to_file(
        bodies, outfile, timestamp=time, append_to_file=False,
    )


def main():
    args = new_argument_parser().parse_args()
    make_oort_cloud(**args.__dict__)


if __name__ == "__main__":
    main()
