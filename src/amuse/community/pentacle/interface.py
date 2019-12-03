"""
Interface for Pentacle
"""

from amuse.rfi.core import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
)
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.community.interface.gd import (
    GravitationalDynamicsInterface,
    GravitationalDynamics,
    GravityFieldInterface,
    GravityFieldCode,
)
from amuse.units import nbody_system


class PentacleInterface(
        CodeInterface,
        LiteratureReferencesMixIn,
        GravitationalDynamicsInterface,
        GravityFieldInterface,
):
    """
    Pentacle hybrid particle-particle particle-tree code based on FDPS

    .. [#] Iwasawa, M. et al. (2016, Publications of the Astronomical Society of Japan, 68, 54)
    .. [#] Namekata, D. et al. (2018, Publications of the Astronomical Society of Japan, 70, 70)
    .. [#] Iwasawa et al., ""
    """
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self,
            name_of_the_worker="pentacle_worker",
            **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def set_time_step():
        """
        Set the current time step.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'time_step', dtype='float64',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_step():
        """
        Get the current time step.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'time_step', dtype='float64',
            direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_eta():
        """
        Set the current accuracy parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'eta', dtype='float64',
            direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        """
        Get the current accuracy parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'eta', dtype='float64',
            direction=function.OUT)
        function.result_type = 'int32'
        return function

class Pentacle(
        GravitationalDynamics,
        GravityFieldCode,
):
    """
    Low-level Pentacle interface
    """

    def __init__(self, convert_nbody=None, **options):
        GravitationalDynamics.__init__(
            self,
            PentacleInterface(**options),
            convert_nbody,
            **options
        )

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_eta",
            "set_eta",
            "time_step_parameter",
            "time step parameter for Hermite integrator",
            default_value = 0.1,
        )
        
        handler.add_method_parameter(
            "get_time_step",
            "set_time_step",
            "time_step",
            "time step for Tree integrator",
            default_value = 0.00390625 | nbody_system.time,
        )

        handler.add_method_parameter(
            "get_eps2",
            "set_eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.01 | nbody_system.length * nbody_system.length,
        )

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            "new_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "set_eps2",
            (
                nbody_system.length * nbody_system.length
            ),
            (
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "get_eps2",
            (),
            (
                nbody_system.length * nbody_system.length,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "set_time_step",
            (
                nbody_system.time,
            ),
            (
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "get_time_step",
            (),
            (
                nbody_system.time,
                handler.ERROR_CODE
            )
        )
