"""
This is the public interface to the AMUSE
*Astrophysical Multipurpose Software Environment* framework.
"""

from amuse.support.core import late

from amuse.datamodel import (
    particle_attributes, Particle, Particles, ParticlesSuperset, Grid,
)
from amuse.io import (
    write_set_to_file, read_set_from_file, get_options_for_format, ReportTable,
)
from amuse.units import (
    units, core, si, constants, nbody_system, generic_unit_system,
    generic_unit_converter,
)
from amuse.units import generic_unit_system as generic
from amuse.units import nbody_system as nbody
from amuse.units.quantities import (
    zero, ScalarQuantity, VectorQuantity, AdaptingVectorQuantity, new_quantity,
)

from amuse.support.console import (
    set_printing_strategy, get_current_printing_strategy,
)

from amuse.ic.plummer import new_plummer_model, new_plummer_sphere
from amuse.ic.salpeter import (
    new_salpeter_mass_distribution, new_salpeter_mass_distribution_nbody,
    new_powerlaw_mass_distribution, new_powerlaw_mass_distribution_nbody,
)
from amuse.ic.brokenimf import (
    new_broken_power_law_mass_distribution, new_scalo_mass_distribution,
    new_miller_scalo_mass_distribution, new_kroupa_mass_distribution,
)
from amuse.ic.flatimf import (
    new_flat_mass_distribution, new_flat_mass_distribution_nbody,
)
from amuse.ic.kingmodel import new_king_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.molecular_cloud import new_ism_cube
from amuse.ext.solarsystem import (
    new_solar_system_for_mercury, new_solar_system,
)
from amuse.ext.spherical_model import (
    new_uniform_spherical_particle_distribution,
    new_spherical_particle_distribution
)


def _placeholder(e):
    class _placeholder(object):
        def __init__(self, *arg, **kwargs):
            raise e
    return _placeholder


try:
    from amuse.ext.halogen_model import new_halogen_model
except ImportError as e:
    locals()["new_halogen_model"] = _placeholder(e)


try:
    from amuse.ext.galactics_model import new_galactics_model
except ImportError as e:
    locals()["new_galactics_model"] = _placeholder(e)


try:
    from amuse.ext.star_to_sph import (
        convert_stellar_model_to_SPH, pickle_stellar_model,
    )
except ImportError as e:
    locals()["convert_stellar_model_to_SPH"] = _placeholder(e)
    locals()["pickle_stellar_model"] = _placeholder(e)


_community_codes = [
    "BHTree",
    "Hermite",
    "PhiGRAPE",
    "Octgrav",
    "TwoBody",
    "Huayno",
    "ph4",
    "Bonsai",
    "Pikachu",
    "AarsethZare",
    "Adaptb",
    "Hacs64",
    "HiGPUs",
    "Kepler",
    "Mercury",
    "MI6",
    "Mikkola",
    "SmallN",
    "Rebound",
    "Brutus",
    "Fi",
    "Gadget2",
    "Athena",
    "Capreole",
    "MpiAmrVac",
    "SimpleX",
    "Mocassin",
    "SPHRay",
    "SSE",
    "BSE",
    "MOSSE",
    "MOBSE",
    "SeBa",
    "EVtwin",
    "Mesa",
    "MMAMS",
    "Hop",
    ]


for _name in _community_codes:
    _interfacename = _name+"Interface"
    # future fix: _interfacename = _name.title()+"Interface"
    _packagename = _name.lower()
    _standardisedname = _name.title()
    try:
        _interface = __import__(
            "amuse.community."+_packagename+".interface",
            fromlist=[_name, _interfacename, _standardisedname],
        )
        locals()[_name] = getattr(_interface, _name)
        locals()[_interfacename] = getattr(_interface, _interfacename)
        locals()[_standardisedname] = getattr(_interface, _standardisedname)
    except ImportError as e:
        locals()[_name] = _placeholder(e)
        locals()[_interfacename] = _placeholder(e)
        locals()[_standardisedname] = _placeholder(e)


def vector(value=[], unit=None):
    if unit is None:
        if isinstance(value, core.unit):
            return VectorQuantity([], unit=value)
        elif isinstance(value, ScalarQuantity):
            return value.as_vector_with_length(1)
        else:
            result = AdaptingVectorQuantity()
            result.extend(value)
            return result
    else:
        if isinstance(value, ScalarQuantity):
            return value.as_vector_with_length(1)
        else:
            return VectorQuantity(value, unit)
