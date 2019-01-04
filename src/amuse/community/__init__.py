import inspect
import itertools
import numpy
import os

from collections import defaultdict

from amuse import datamodel
from amuse.community.interface.stopping_conditions import (
    StoppingConditionInterface, StoppingConditions,
)
from amuse.datamodel import (
    base, parameters, incode_storage, attributes,
)
from amuse.rfi.core import (
    CodeInterface, CodeWithDataDirectories, legacy_function, remote_function,
    LegacyFunctionSpecification, is_mpd_running,
)
from amuse.support import (
    exceptions, state, get_amuse_root_dir,
)
from amuse.support.core import late
from amuse.support.interface import (
    ConvertArgumentsException, OldObjectsBindingMixin,
    MethodArgumentOrResultType, NoUnitMethodArgumentOrResultType,
    UnitMethodArgumentOrResultType, ErrorCodeMethodArgumentOrResultType,
    IndexMethodArgumentOrResultType, _get_result_type,
    LinkMethodArgumentOrResultType, CodeAttributeWrapper,
    HandleCodeInterfaceAttributeAccess, LegacyInterfaceHandler,
    HandleConvertUnits, StateMethodDefinition, HandleState, MethodWithUnits,
    MethodWithUnitsDefinition, HandleMethodsWithUnits,
    PropertyWithUnitsDefinition, HandlePropertiesWithUnits, HandleParameters,
    HandleErrorCodes, AbstractParticleSetDefinition, ParticleSetDefinition,
    ParticleSupersetDefinition, GridDefinition, CodeInMemoryParticles,
    HandleParticles, OverriddenCodeInterface, InCodeComponentImplementation,
    IncorrectMethodDefinition, PropertyDefinition,
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.support.methods import (
    CodeMethodWrapper, CodeMethodWrapperDefinition,
    IncorrectWrappedMethodException, ProxyingMethodWrapper,
)
from amuse.support.options import OptionalAttributes, option
from amuse.units import (
    nbody_system, generic_unit_system, quantities, core, units,
)
from amuse.units.core import unit
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.units.quantities import is_quantity


ERROR_CODE = MethodWithUnitsDefinition.ERROR_CODE
NO_UNIT = MethodWithUnitsDefinition.NO_UNIT
INDEX = MethodWithUnitsDefinition.INDEX
LINK = MethodWithUnitsDefinition.LINK

"""
Existing, production codes

Contains the source code of production codes and software to embed these codes
into AMUSE
"""
