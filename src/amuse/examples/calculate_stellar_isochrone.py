"""
   Evolve a population of N stars.
   Initial mass function between Mmin and Mmax with stellar evolution
   for metallicity z.
"""
import numpy as np
from amuse.io import write_set_to_file, read_set_from_file
from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.seba import Seba
from amuse.community.sse import Sse
from amuse.community.mesa_r2208.interface import Mesa
#from amuse.community.mesa import Mesa
from amuse.community.evtwin import Evtwin
from amuse.ic.salpeter import new_salpeter_mass_distribution

def calculate_stellar_temperature_and_luminosity(
        star,
        evo_code, metallicity=0.02, time_end=100 | units.Myr):

    if evo_code.find("Seba") >= 0:
        stellar = Seba()
        stellar.parameters.metallicity = metallicity
        stellar.particles.add_particle(star)
        channel_to_framework = stellar.particles.new_channel_to(star)
        stellar.evolve_model(time_end)
        channel_to_framework.copy_attributes(
            ["radius", "temperature", "luminosity"]
        )
        stellar.stop()
    else:
        if evo_code.find("Mesa") >= 0:
            stellar = Mesa()
        elif evo_code.find("Evtwin") >= 0:
            stellar = Evtwin()
        stellar.parameters.metallicity = metallicity
        stellar.particles.add_particle(star)
        channel_to_framework = stellar.particles.new_channel_to(star)
        try:
            if star.mass<=0.8|units.MSun:
                t_end = 0.1|units.Myr
            else:
                t_end = time_end
            print("Evolved star: m=", star.mass.in_(units.MSun))
            stellar.evolve_model(t_end)
            channel_to_framework.copy_attributes(
                ["radius", "temperature", "luminosity"])
            print("Successfully evolved star: m=", star.mass.in_(units.MSun))
        except:
            print("Failed to evolve star: m=", star.mass.in_(units.MSun))
        stellar.stop()
    return star

def main(input_file, time_end, metallicity, evo_code="Seba"):
    evo_code = evo_code.capitalize()
    if evo_code.find("Seba") >= 0:
        filename = "../data/Stellar_Isochrone_Seba.amuse"
    elif evo_code.find("Mesa") >= 0:
        filename = "../data/Stellar_Isochrone_Mesa.amuse"
    elif evo_code.find("Evtwin") >= 0:
        filename = "../data/Stellar_Isochrone_Evtwin.amuse"
    else:
        print("No input code, stop")
        exit(-1)
        
    stars = read_set_from_file(input_file, "amuse", close_file=True)
    from os.path import exists
    if exists(filename):
        stars_done = read_set_from_file(filename, "amuse", close_file=True)
    else:
        stars_done = Particles(0)
    for si in stars:
        if si not in stars_done:
            star_evolved = calculate_stellar_temperature_and_luminosity(
                si.as_set(),
                evo_code=evo_code, metallicity=metallicity,
                time_end=time_end)
            stars_done.add_particle(star_evolved)
            write_set_to_file(stars_done, filename, overwrite_file=True)


def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-C", dest="evo_code", default="Seba",
                      help="stellar evolution code [SeBa]")
    result.add_option("--f", dest="input_file",
                      default="initial_stellar_mass_function.amuse",
                      help="input data file")
    result.add_option("-t", dest="time_end", unit=units.Myr,
                      type="float", default=4500.0 | units.Myr,
                      help="end time of the simulation [4500] Myr")
    result.add_option("-z", dest="metallicity", type="float", default=0.02,
                      help="metalicity [0.02]")
    return result


if __name__ in ('__main__'):
    o, arguments = new_option_parser().parse_args()
    main(**o.__dict__)
