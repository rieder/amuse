import warnings


def rename_fn_par(new_name, new_value, old_name, old_value, default_value):
    """Get value of a renamed function parameter.

    If you have a function f with parameter ``x`` and default value 1, like this:

    .. code-block::

        def f(x=1):
            ...

    and you want to rename ``x`` to ``y`` while staying backwards compatible, then you
    can rename ``x`` to ``y``, add ``x`` back in at the end so that old keyword
    arguments still work, give both a default value of ``None``, and use this function
    like this:

    .. code-block::

        def f(y=None, x=None):
            value = rename_fn_par("y", y, "x", x, 1)

    Any callers using ``x`` explicitly will receive a warning to change their code to
    use ``y`` in the future. If both ``x`` and ``y`` are set and the values are
    different, then an exception is raised.

    Args:
        new_name (str): The new name of the variable
        new_value: The value passed using the new name
        old_name (str): The old name of the variable
        old_value: The value passed using the old name
        default_value: The default value if neither are set

    Returns:
        Either new_value or old_value if only one is set or they're set to the same
        value, or default_value if neither is set.

    Raises:
        ValueError: If both new_value and old_value are set, and to different values.
    """
    if new_value is not None:
        if old_value is not None and old_value != new_value:
            raise ValueError(
                    f"{old_name} and {new_name} have different values,"
                    " which is not allowed because they represent the same thing.")
        return new_value

    if old_value is not None:
        warnings.warn(
                f"{old_name} is deprecated, please use {new_name} instead",
                category=FutureWarning)
        return old_value

    return default_value
