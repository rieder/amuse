"""
Runs AMUSE examples (or a selection thereof), and benchmarks them.
"""
import sys
import time
import os
import argparse
import platform
import contextlib
if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files
from datetime import datetime, timezone
import matplotlib.pyplot as plt
from amuse.units import units
import amuse.examples as ae
from .list_of_examples import examples_per_chapter

#  uname_result(system='Linux', node='braam', release='6.12.20+rpt-rpi-2712', version='#1 SMP PREEMPT Debian 1:6.12.20-1+rpt1~bpo12+1 (2025-03-19)', machine='aarch64')
def run_example(name, restart=True):
    if restart:
        if os.path.isfile(f"{name}.out"):
            os.remove(f"{name}.out")
        if os.path.isfile(f"{name}.time"):
            os.remove(f"{name}.time")
    else:
        if os.path.isfile(f"{name}.out"):
            print(f"Skipping example {name}")
            return
    start_time = time.time() | units.s
    with open(f"{name}.out", "w") as f:
        with contextlib.redirect_stdout(f):
            try:
                ae.run(name)
                success = True
            except Exception as ex:
                sys.stderr.write(f"Failed to run example {name}\n")
                sys.stderr.write(f"{ex}\n")
                success = False

    if not success:
        os.remove(f"{name}.out")
        return
    end_time = time.time() | units.s
    duration = end_time - start_time

    timefile = f"{name}.time"
    if not os.path.isfile(timefile):
        with open(timefile, "w", encoding="utf-8") as out:
            out.write("# duration\ttime of execution\tplatform used\tarchitecture used\tPython version used\n")
    with open(f"{name}.time", "a", encoding="utf-8") as out:
        out.write(f"{duration.in_(units.s)}\t{datetime.now(timezone.utc)}\t{platform.system()}\t{platform.machine()}\t{platform.python_version()}\n")
    return


def new_argument_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-e", "--example",
        nargs='+',
        default=None,
        help="Name of the example(s) to run",
    )

    parser.add_argument(
        "-c", "--chapter",
        nargs='+',
        default=None,
        help="Chapter(s) of which all examples should be run",
    )

    parser.add_argument(
        "-r", "--restart",
        action="store_true",
        default=False,
        help="Restart benchmarking",
    )
    parser.add_argument(
        "-p",
        "--plot_style",
        type=str,
        default="amusebook.mplstyle",
        help="Matplotlib style file to use",
    )
    return parser


def main():
    arguments = new_argument_parser().parse_args()
    plotstylefile = files("amuse.examples.data").joinpath(arguments.plot_style)
    plt.style.use(plotstylefile)

    examples_to_run = set()
    print(arguments.chapter)
    print(arguments.example)
    if arguments.chapter is not None:
        for ch in arguments.chapter:
            if ch in examples_per_chapter:
                for ex in examples_per_chapter[ch]:
                    examples_to_run.add(ex)
    if arguments.example is not None:
        for ex in arguments.example:
            examples_to_run.add(ex)

    for ex in examples_to_run:
        print(f"Running example {ex}")
        try:
            run_example(ex, arguments.restart)
        except Exception as fail:
            print(f"Failed to run example {ex}: {fail}")
    return


if __name__ == "__main__":
    main()
