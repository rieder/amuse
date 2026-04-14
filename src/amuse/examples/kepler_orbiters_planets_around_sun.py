"""
simple example to use kepler_orbiters
integrates the solar system planets motion around the Sun
"""

from amuse.units import units, nbody_system, quantities
from amuse.datamodel import Particles
from amuse.community.kepler_orbiters.interface import Kepler
from matplotlib import pyplot


def sun_and_planets():
    """
    solar system on 2014-Jul-09 (JD 2456847.5)
      Sun + 8 planets
      http://ssd.jpl.nasa.gov/horizons.cgi
    """
    sun = Particles(1)
    sun[0].mass = 1.0 | units.MSun
    sun[0].position = (
        2.037427757106262e-03,
        -1.670788401563819e-03,
        -1.177108759951928e-04,
    ) | units.au
    sun[0].velocity = (
        5.088355581933507e-06,
        3.932938906525153e-06,
        -1.245782562840000e-07,
    ) | (units.au / units.day)
    sun[0].radius = 1.0 | units.RSun

    planets = Particles(8)
    mercury = planets[0]
    mercury.mass = 3.302e23 | units.kg
    mercury.position = (
        3.387829599522694e-01,
        -2.068479571871339e-01,
        -4.777798825282107e-02,
    ) | units.au
    mercury.velocity = (
        9.121685164224357e-03,
        2.532046488259770e-02,
        1.231999349286273e-03,
    ) | (units.au / units.day)
    mercury.radius = 2440.0 | units.km

    venus = planets[1]
    venus.mass = 48.685e23 | units.kg
    venus.position = (
        5.829554745446001e-01,
        4.290508478691337e-01,
        -2.773923286024579e-02,
    ) | units.au
    venus.velocity = (
        -1.210377003565980e-02,
        1.616561583873175e-02,
        9.201724371252968e-04,
    ) | (units.au / units.day)
    venus.radius = 6051.8 | units.km

    earth = planets[2]
    earth.mass = 5.97219e24 | units.kg
    earth.position = (
        2.914170806309638e-01,
        -9.762380901507695e-01,
        -8.838036341545362e-05,
    ) | units.au
    earth.velocity = (
        1.621194321564609e-02,
        4.839618876688764e-03,
        -8.288236310894908e-07,
    ) | (units.au / units.day)
    earth.radius = 6371.01 | units.km

    mars = planets[3]
    mars.mass = 6.4185e23 | units.kg
    mars.position = (
        -6.771052080404090e-01,
        -1.358079217192958e00,
        -1.186843692826959e-02,
    ) | units.au
    mars.velocity = (
        1.304721335990521e-02,
        -5.060040804946996e-03,
        -4.263506287278560e-04,
    ) | (units.au / units.day)
    mars.radius = 3389.9 | units.km

    jupiter = planets[4]
    jupiter.mass = 1898.13e24 | units.kg
    jupiter.position = (
        -2.659291255441077e00,
        4.536415692760404e00,
        4.058673324381942e-02,
    ) | units.au
    jupiter.velocity = (
        -6.601223938431606e-03,
        -3.460147408421039e-03,
        1.620752232688015e-04,
    ) | (units.au / units.day)
    jupiter.radius = 71492.0 | units.km

    saturn = planets[5]
    saturn.mass = 5.68319e26 | units.kg
    saturn.position = (
        -6.149519345675703e00,
        -7.775694595902171e00,
        3.799392627119177e-01,
    ) | units.au
    saturn.velocity = (
        4.071035079488449e-03,
        -3.476237592282034e-03,
        -1.012547571724344e-04,
    ) | (units.au / units.day)
    saturn.radius = 60268.0 | units.km

    uranus = planets[6]
    uranus.mass = 86.8103e24 | units.kg
    uranus.position = (
        1.948219896924911e01,
        4.611924353502811e00,
        -2.352696714787791e-01,
    ) | units.au
    uranus.velocity = (
        -9.347014358359851e-04,
        3.643979378812646e-03,
        2.560171100139048e-05,
    ) | (units.au / units.day)
    uranus.radius = 25559.0 | units.km

    neptune = planets[7]
    neptune.mass = 102.41e24 | units.kg
    neptune.position = (
        2.731034114396837e01,
        -1.235250519790972e01,
        -3.750181405673675e-01,
    ) | units.au
    neptune.velocity = (
        1.272538140363074e-03,
        2.879023382710963e-03,
        -8.821000085244238e-05,
    ) | (units.au / units.day)
    neptune.radius = 24766.0 | units.km

    # coordinates transformation
    # with respect to solar system barycenter --> with respect to the Sun
    sun.position -= sun.position
    sun.velocity -= sun.velocity
    planets.position -= sun.position
    planets.velocity -= sun.velocity

    return sun, planets


def integrate_solar_system(sun, planets, time_end=5.0 | units.yr, n_steps=500):
    """
    evolve the system using kepler_orbiters
    """

    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.au)

    planets_around_sun = Kepler(converter, channel_type="sockets")

    # central particle
    planets_around_sun.central_particle.add_particles(sun[0:1])

    # to set the central particle at the center of the coordinate system
    # planets_around_sun.central_particle.position = (0.0, 0.0, 0.0) | units.au
    # planets_around_sun.central_particle.velocity = (0.0, 0.0, 0.0) | units.kms

    # orbiters
    planets_around_sun.orbiters.add_particles(planets)
    planets_around_sun.commit_particles()

    # to change the integration method
    # planets_around_sun.parameters.method = 1
    # print planets_around_sun.get_method()

    channel_from_planetets_to_framework = planets_around_sun.orbiters.new_channel_to(
        planets
    )
    channel_from_sun_to_framework = planets_around_sun.central_particle.new_channel_to(
        sun
    )

    positions_sun = quantities.AdaptingVectorQuantity()
    positions_planets = quantities.AdaptingVectorQuantity()

    dt = time_end / float(n_steps)
    time = 0.0 | units.yr

    print(
        " ** evolving solar system for",
        time_end.in_(units.yr),
        ", with time-step of",
        dt.in_(units.yr),
    )
    print("    this might take a while")
    while time <= time_end:
        # print "\t", time.in_(units.yr)
        planets_around_sun.evolve_model(time)
        channel_from_planetets_to_framework.copy()
        channel_from_sun_to_framework.copy()
        positions_sun.append(sun.position)
        positions_planets.append(planets.position)
        time += dt
    print(" **")

    planets_around_sun.stop()

    return positions_sun, positions_planets


def plot_solar_system(positions_sun, positions_planets):
    """
    simple plot to get the XY and XZ plane
    """

    fig = pyplot.figure(figsize=(10, 7))

    # XY plane
    pyplot.subplot(121)
    pyplot.scatter(
        positions_sun[0, 0, 0].value_in(units.au),
        positions_sun[0, 0, 1].value_in(units.au),
    )
    for i, planet_i in enumerate(positions_planets[0, :, 0]):
        pyplot.plot(
            positions_planets[:, i, 0].value_in(units.au),
            positions_planets[:, i, 1].value_in(units.au),
        )
    pyplot.axis([-45, 45, -45, 45])
    pyplot.gca().set_aspect("equal")
    pyplot.xlabel("x [AU]")
    pyplot.ylabel("y [AU]")

    # XZ plane
    pyplot.subplot(122)
    pyplot.scatter(
        positions_sun[0, 0, 0].value_in(units.au),
        positions_sun[0, 0, 2].value_in(units.au),
    )
    for i, planet_i in enumerate(positions_planets[0, :, 0]):
        pyplot.plot(
            positions_planets[:, i, 0].value_in(units.au),
            positions_planets[:, i, 2].value_in(units.au),
        )
    pyplot.axis([-45, 45, -10, 10])
    pyplot.gca().set_aspect("equal")
    pyplot.xlabel("x [au]")
    pyplot.ylabel("z [au]")

    pyplot.tight_layout()
    pyplot.show()

    return


if __name__ == "__main__":

    # initial conditions
    sun, planets = sun_and_planets()

    # integrate the solar system
    positions_sun, positions_planets = integrate_solar_system(
        sun, planets, time_end=10.0 | units.yr, n_steps=1000
    )

    # simple plot of the results
    plot_solar_system(positions_sun, positions_planets)
