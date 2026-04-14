from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.seba import Seba
from amuse.community.mesa import Mesa
from amuse.community.sse import Sse
from amuse.support.console import set_printing_strategy


# #BOOKLISTSTART# #
def evolve_with_different_stellar_model():
    times = [10, 100, 1000] | units.Myr
    stars = Particles(mass=[1, 2, 4] | units.MSun)
    stellars = [Seba(), Mesa(), Sse()]
    channels = []
    for star, stellar in zip(stars, stellars):
        stellar.particles.add_particle(star)
        channels.append(stellar.particles.new_channel_to(stars))

    for time, channel, stellar in zip(times, channels, stellars):
        stellar.evolve_model(time)
        channel.copy()

    for time, star in zip(times, stars):
        print(f"Time= {time} M= {star.mass}")

    for stellar in stellars:
        stellar.stop()


# #BOOKLISTSTOP# #


def evolve_with_same_stellar_model():
    time = [10, 100, 1000] | units.Myr
    stars = Particles(mass=[1, 2, 4] | units.MSun)
    stellar = Seba()
    stellar.particles.add_particles(stars)
    channel = stellar.particles.new_channel_to(stars)
    for ti, si in zip(time, stellar.particles):
        si.evolve_for(ti)
        channel.copy()

    print(f"Time= {stellar.model_time} age= {stars.age}\n" f"\tmass= {stars.mass}")
    stellar.stop()


def multiple_stellar_code():
    evolve_with_different_stellar_model()
    evolve_with_same_stellar_model()


def main(**kwargs):
    multiple_stellar_code(**kwargs)


if __name__ in ("__main__", "__plot__"):
    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=6,
        prefix="",
        separator="[",
        suffix="]",
    )
    main()
