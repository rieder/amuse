from amuse.community import (
        CodeInterface,
        LegacyFunctionSpecification,
        legacy_function,
        LiteratureReferencesMixIn,
        # StoppingConditionInterface,
        )
from amuse.community.interface.gd import (
        GravitationalDynamicsInterface,
        GravitationalDynamics,
        GravityFieldInterface,
        GravityFieldCode,
        )


class GasolineInterface(
        CodeInterface,
        GravitationalDynamicsInterface,
        LiteratureReferencesMixIn,
        # StoppingConditionInterface,
        # GravityFieldInterface,
        # SinglePointGravityFieldInterface,
        ):
    """
    Gasoline: a flexible, parallel implementation of TreeSPH
    

    The relevant references are:
        .. [#] Wadsley, Keller, Quinn, 2017, MNRAS 471, 2357
        .. [#] Wadsley, Stadel, Quinn, 2004, NewA 9, 137
    """

    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
                self,
                name_of_the_worker="gasoline_worker",
                **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    def new_particle(self, mass, x, y, z, vx, vy, vz, u=0.):
        return self.new_sph_particle(mass, x, y, z, vx, vy, vz, u)

    @legacy_function
    def new_sph_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.OUT)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x', 'y', 'z']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x', 'y', 'z']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def set_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter(
                'number_of_particles',
                dtype='i',
                direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function


class Gasoline(
        GravitationalDynamics,
        GravityFieldCode,
        ):

    def __init__(self, convert_nbody=None, **options):

        # self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
                self,
                GasolineInterface(**options),
                convert_nbody,
                **options)
