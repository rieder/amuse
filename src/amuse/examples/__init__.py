"""
AMUSE examples module
"""

import inspect
import importlib

try:
    from IPython.core.getipython import get_ipython
    from IPython.terminal.interactiveshell import InteractiveShell
except ImportError:

    def get_ipython():
        "N/A, so return None"
        return None


def get(example_name, prefix="amuse_examples"):
    """
    Retrieves an AMUSE example.
    """
    if prefix is not None:
        if prefix[-1] != ".":
            prefix += "."
    return inspect.getsource(importlib.import_module(f"{prefix}{example_name}"))


def to_cell(example_name, prefix="amuse_examples"):
    """
    Creates a Jupyter cell from an AMUSE example.
    Works also for other modules, set prefix to None in that case.
    """

    contents = get(example_name, prefix=prefix)
    shell = get_ipython()
    shell.set_next_input(contents, replace=False)


def show(example_name, prefix="amuse_examples"):
    """
    Prints the source code of an AMUSE example.
    """
    print(get(example_name, prefix=prefix))


def run(example_name):
    """
    Runs an AMUSE example.
    """
    prefix = "amuse_examples."
    try:
        example = importlib.import_module(f"{prefix}{example_name}")
    except ImportError as exc:
        raise ImportError("Could not find example") from exc
    example.main()


def shell_is_interactive():
    """
    Returns True if the shell is interactive and False otherwise.
    """
    try:
        return isinstance(get_ipython(), InteractiveShell)
    except ImportError:
        return False
