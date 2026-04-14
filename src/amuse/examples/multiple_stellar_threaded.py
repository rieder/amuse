"""
Evolve stars with multiple stellar codes, threaded version
"""

import sys
import argparse
import queue
import threading
import multiprocessing

from amuse.datamodel import Particles
from amuse.units import units
from amuse.support.console import set_printing_strategy

try:
    from amuse.community.sse import Sse
except ModuleNotFoundError:
    Sse = None
try:
    from amuse.community.seba import Seba
except ModuleNotFoundError:
    Seba = None
try:
    from amuse.community.evtwin import Evtwin
except ModuleNotFoundError:
    Evtwin = None
try:
    from amuse.community.mesa import Mesa
except ModuleNotFoundError:
    Mesa = None

# #BOOKLISTSTART1# #
code_queue = queue.Queue()


def remote_worker_code():
    code = code_queue.get()
    evolve_single_star(code)
    code_queue.task_done()


def evolve_with_different_stellar_model(codes):
    for code in codes:
        if code is not None:
            code_queue.put(code)
    n_cpu = multiprocessing.cpu_count()
    for i in range(n_cpu):
        thread = threading.Thread(target=remote_worker_code)
        thread.daemon = True
        thread.start()
    code_queue.join()  # block until all tasks are done


# #BOOKLISTSTOP1# #


# #BOOKLISTSTART2# #
def evolve_single_star(code):
    stars = Particles(mass=10 | units.MSun)
    stellar = code()
    stellar.particles.add_particles(stars)
    channel = stellar.particles.new_channel_to(stars)

    stellar.evolve_model(1 | units.Myr)
    channel.copy()
    print(
        f"Star evolved to time= {stellar.model_time} "
        f"M= {stars.mass} R= {stars.radius}"
    )
    stellar.stop()


# #BOOKLISTSTOP2# #


def new_argument_parser():
    result = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument("-t", action="store_true", dest="threaded", help="run threaded")
    return result


def main(args=sys.argv[1:]):
    arguments = new_argument_parser().parse_args(args)
    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.RSun, units.Myr],
        precision=6,
        prefix="",
        separator="[",
        suffix="]",
    )

    codes = [Seba, Mesa, Sse, Evtwin]

    if arguments.threaded:
        print("Run threaded")
        evolve_with_different_stellar_model(codes)
    else:
        print("Run sequentially")
        for code in codes:
            if code is not None:
                evolve_single_star(code)


if __name__ == "__main__":
    main()
