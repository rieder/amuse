"Interface to Phantom"

from amuse.community import (
    CodeInterface,
    LegacyFunctionSpecification,
    legacy_function,
    LiteratureReferencesMixIn,
)

from amuse.community.interface.gd import (
    GravitationalDynamicsInterface,
    GravitationalDynamics,
    # GravityFieldInterface,
    GravityFieldCode,
)
from amuse.community.interface.stopping_conditions import(
    StoppingConditionInterface,
    StoppingConditions,
)
from amuse.units import nbody_system, units


class PhantomInterface(
        CodeInterface,
        LiteratureReferencesMixIn,
        GravitationalDynamicsInterface,
        StoppingConditionInterface,
        # SinglePointGravityFieldInterface,
):
    """
    The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al.

    References:
        .. [#] Price et al., 2018, PASA, Volume 35, id.e031 82 pp
    """

    def __init__(self, **options):
        CodeInterface.__init__(
            self,
            name_of_the_worker="phantom_worker",
            **options)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def new_dm_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
        )
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter(
            'radius', dtype='float64', direction=function.IN, default=0.01,
        )
        function.result_type = 'int32'
        return function

    def new_particle(self, mass, x, y, z, vx, vy, vz, radius):
        return self.new_dm_particle(mass, x, y, z, vx, vy, vz, radius)

    @legacy_function
    def new_sph_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
        )
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN, default=0.01,
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_sink_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.OUT,
        )
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter(
            'radius', dtype='float64', direction=function.IN, default=0.01,
            # default should be h_acc
        )
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN, default=0.01,
            # default should be h_smooth_sinksink?
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state_dm():
        """
        Retrieve the current state of a DM particle. The mass, position and
        velocity are returned.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle to get the state from. This
            index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.OUT,
            description="The current mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'radius', dtype='float64', direction=function.OUT,
            description="The current softening length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_state_sink():
        """
        Retrieve the current state of a sink particle. The mass, position and
        velocity are returned.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle to get the state from. This
            index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.OUT,
            description="The current mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'radius', dtype='float64', direction=function.OUT,
            description="The accretion radius of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="The smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    def get_state(self, index_of_the_particle):
        return self.get_state_dm(index_of_the_particle)

    @legacy_function
    def get_state_sph():
        """
        Retrieve the current state of an SPH particle. The mass, position,
        velocity, internal energy and smoothing length are returned.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle to get the state from. This
            index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.OUT,
            description="The current mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.OUT,
            description="The current position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.OUT,
            description="The current velocity vector of the particle")
        function.addParameter(
            'u', dtype='float64', direction=function.OUT,
            description="The current internal energy of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="The current smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_state_sph():
        """
        Update the current state of an SPH particle. The mass, position,
        velocity, internal energy and smoothing length are updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle for which the state is to be
            updated. This index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The new mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'u', dtype='float64', direction=function.IN,
            description="The new internal energy of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description="The new smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    def set_state(self, index_of_the_particle, mass, x, y, z, vx, vy, vz, radius):
        return self.set_state_dm(
            index_of_the_particle, mass, x, y, z, vx, vy, vz, radius)

    @legacy_function
    def set_state_dm():
        """
        Update the current state of a DM particle. The mass, position and
        velocity are updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle for which the state is to be
            updated. This index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The new mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'radius', dtype='float64', direction=function.IN,
            description="The new softening length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def set_state_sink():
        """
        Update the current state of a sink particle. The mass, position and
        velocity are updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description="""Index of the particle for which the state is to be
            updated. This index must have been returned by an earlier call to
            :meth:`new_particle`""")
        function.addParameter(
            'mass', dtype='float64', direction=function.IN,
            description="The new mass of the particle")
        function.addParameter(
            'x', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'y', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'z', dtype='float64', direction=function.IN,
            description="The new position vector of the particle")
        function.addParameter(
            'vx', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vy', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'vz', dtype='float64', direction=function.IN,
            description="The new velocity vector of the particle")
        function.addParameter(
            'radius', dtype='float64', direction=function.IN,
            description="The accretion radius of the particle")
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description="The smoothing length of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def set_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description='',
        )
        function.addParameter(
            'u', dtype='float64', direction=function.IN,
            description='',
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_smoothing_length():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description='',
        )
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.IN,
            description=''
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=''
        )
        function.addParameter(
            'density', dtype='float64', direction=function.OUT,
            description="The current density of the particle"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_pressure():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description='',
        )
        function.addParameter(
            'pressure', dtype='float64', direction=function.OUT,
            description=''
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=''
        )
        function.addParameter(
            'u', dtype='float64', direction=function.OUT,
            description="The current internal energy of the particle",
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_smoothing_length():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=''
        )
        function.addParameter(
            'h_smooth', dtype='float64', direction=function.OUT,
            description="The current smoothing length of the particle"
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_c_courant():
        function = LegacyFunctionSpecification()
        function.addParameter('c_courant', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_c_courant():
        function = LegacyFunctionSpecification()
        function.addParameter('c_courant', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_dtmax():
        function = LegacyFunctionSpecification()
        function.addParameter('dtmax', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_dtmax():
        function = LegacyFunctionSpecification()
        function.addParameter('dtmax', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_c_force():
        function = LegacyFunctionSpecification()
        function.addParameter('C_force', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_c_force():
        function = LegacyFunctionSpecification()
        function.addParameter('C_force', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_tolv():
        function = LegacyFunctionSpecification()
        function.addParameter('tolv', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_tolv():
        function = LegacyFunctionSpecification()
        function.addParameter('tolv', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_hfact():
        function = LegacyFunctionSpecification()
        function.addParameter('hfact', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_hfact():
        function = LegacyFunctionSpecification()
        function.addParameter('hfact', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_tolh():
        function = LegacyFunctionSpecification()
        function.addParameter('tolh', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_tolh():
        function = LegacyFunctionSpecification()
        function.addParameter('tolh', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_tree_accuracy():
        function = LegacyFunctionSpecification()
        function.addParameter('tree_accuracy', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_tree_accuracy():
        function = LegacyFunctionSpecification()
        function.addParameter('tree_accuracy', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_alpha():
        function = LegacyFunctionSpecification()
        function.addParameter('alpha', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_alpha():
        function = LegacyFunctionSpecification()
        function.addParameter('alpha', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_alphamax():
        function = LegacyFunctionSpecification()
        function.addParameter('alphamax', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_alphamax():
        function = LegacyFunctionSpecification()
        function.addParameter('alphamax', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_beta():
        function = LegacyFunctionSpecification()
        function.addParameter('beta', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_beta():
        function = LegacyFunctionSpecification()
        function.addParameter('beta', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_avdecayconst():
        function = LegacyFunctionSpecification()
        function.addParameter('avdecayconst', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_avdecayconst():
        function = LegacyFunctionSpecification()
        function.addParameter('avdecayconst', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_idamp():
        function = LegacyFunctionSpecification()
        function.addParameter('idamp', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_idamp():
        function = LegacyFunctionSpecification()
        function.addParameter('idamp', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_ieos():
        function = LegacyFunctionSpecification()
        function.addParameter('ieos', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_ieos():
        function = LegacyFunctionSpecification()
        function.addParameter('ieos', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_gamma():
        function = LegacyFunctionSpecification()
        function.addParameter('gamma', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_mu():
        function = LegacyFunctionSpecification()
        function.addParameter('mu', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_mu():
        function = LegacyFunctionSpecification()
        function.addParameter('mu', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_rhofinal():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'rhofinal', dtype='float64', direction=function.OUT,
            unit=(nbody_system.density),
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_rhofinal():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'rhofinal', dtype='float64', direction=function.IN,
            unit=(nbody_system.density),
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_rho_crit():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'rho_crit', dtype='float64', direction=function.OUT,
            unit=(nbody_system.density),
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_rho_crit():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'rho_crit', dtype='float64', direction=function.IN,
            unit=(nbody_system.density),
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_r_crit():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_crit', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_r_crit():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'r_crit', dtype='float64', direction=function.IN,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_h_acc():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'h_acc', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_h_acc():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'h_acc', dtype='float64', direction=function.IN,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_h_soft_sinkgas():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'h_soft_sinkgas', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_h_soft_sinkgas():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'h_soft_sinkgas', dtype='float64', direction=function.IN,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_h_soft_sinksink():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'h_soft_sinksink', dtype='float64', direction=function.OUT,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_h_soft_sinksink():
        function = LegacyFunctionSpecification()
        function.addParameter(
            'h_soft_sinksink', dtype='float64', direction=function.IN,
            unit=nbody_system.length,
        )
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_f_acc():
        function = LegacyFunctionSpecification()
        function.addParameter('f_acc', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_f_acc():
        function = LegacyFunctionSpecification()
        function.addParameter('f_acc', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_iexternalforce():
        function = LegacyFunctionSpecification()
        function.addParameter('iexternalforce', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_iexternalforce():
        function = LegacyFunctionSpecification()
        function.addParameter('iexternalforce', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_irealvisc():
        function = LegacyFunctionSpecification()
        function.addParameter('irealvisc', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_irealvisc():
        function = LegacyFunctionSpecification()
        function.addParameter('irealvisc', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_shearparam():
        function = LegacyFunctionSpecification()
        function.addParameter('shearparam', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_shearparam():
        function = LegacyFunctionSpecification()
        function.addParameter('shearparam', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def get_bulkvisc():
        function = LegacyFunctionSpecification()
        function.addParameter('bulkvisc', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    @legacy_function
    def set_bulkvisc():
        function = LegacyFunctionSpecification()
        function.addParameter('bulkvisc', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
        -1 - ERROR
        """
        return function

    def reinitialize_particles(self):
        self.recommit_particles()


class Phantom(GravitationalDynamics, GravityFieldCode):
    __interface__ = PhantomInterface

    def __init__(
            self,
            convert_nbody=None,
            **options):
        self.stopping_conditions = StoppingConditions(self)

        GravitationalDynamics.__init__(
            self,
            PhantomInterface(**options),
            convert_nbody,
            **options
        )

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)
        GravityFieldCode.define_state(self, handler)
        # self.stopping_conditions.define_state(handler)

        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)
        handler.add_method('END', 'initialize_code')

        handler.add_transition('RUN', 'UPDATE', 'new_sph_particle', False)
        handler.add_method('EDIT', 'new_sph_particle')
        handler.add_method('UPDATE', 'new_sph_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_dm_particle', False)
        handler.add_method('EDIT', 'new_dm_particle')
        handler.add_method('UPDATE', 'new_dm_particle')
        handler.add_transition('RUN', 'UPDATE', 'new_sink_particle', False)
        handler.add_method('EDIT', 'new_sink_particle')
        handler.add_method('UPDATE', 'new_sink_particle')

        self.stopping_conditions.define_state(handler)

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_c_courant",
            "set_c_courant",
            "C_cour",
            "Courant number",
            default_value=0.3
        )

        handler.add_method_parameter(
            "get_c_force",
            "set_c_force",
            "C_force",
            "dt_force number",
            default_value=0.25
        )

        handler.add_method_parameter(
            "get_tolv",
            "set_tolv",
            "tolv",
            "tolerance on v iterations in timestepping",
            default_value=1.0e-2
        )

        handler.add_method_parameter(
            "get_hfact",
            "set_hfact",
            "hfact",
            "h in units of particle spacing [h = hfact(m/rho)^(1/3)]",
            default_value=1.2
        )

        handler.add_method_parameter(
            "get_tolh",
            "set_tolh",
            "tolh",
            "tolerance on h-rho iterations",
            default_value=1.0e-4
        )

        handler.add_method_parameter(
            "get_tree_accuracy",
            "set_tree_accuracy",
            "tree_accuracy",
            "tree opening criterion (0.0-1.0)",
            default_value=0.5
        )

        handler.add_method_parameter(
            "get_alpha",
            "set_alpha",
            "alpha",
            "MINIMUM art. viscosity parameter",
            default_value=0.
        )

        handler.add_method_parameter(
            "get_alphamax",
            "set_alphamax",
            "alphamax",
            "MAXIMUM art. viscosity parameter",
            default_value=1.0
        )

        handler.add_method_parameter(
            "get_beta",
            "set_beta",
            "beta",
            "beta viscosity",
            default_value=2.0
        )

        handler.add_method_parameter(
            "get_avdecayconst",
            "set_avdecayconst",
            "avdecayconst",
            "decay time constant for viscosity switches",
            default_value=0.1
        )

        handler.add_method_parameter(
            "get_idamp",
            "set_idamp",
            "idamp",
            "artificial damping of velocities (0=off, 1=constant, 2=star)",
            default_value=0
        )

        handler.add_method_parameter(
            "get_ieos",
            "set_ieos",
            "ieos",
            "eqn of state (1=isoth;2=adiab;3=locally iso;8=barotropic)",
            default_value=1
        )

        handler.add_method_parameter(
            "get_gamma",
            "set_gamma",
            "gamma",
            "gamma value ",
            default_value=1
        )

        handler.add_method_parameter(
            "get_mu",
            "set_mu",
            "mu",
            "mean molecular weight",
            default_value=2.381
        )

        handler.add_method_parameter(
            "get_rhofinal",
            "set_rhofinal",
            "rhofinal",
            "maximum allowed density (<=0 to ignore)",
            default_value=(0 | nbody_system.density)
        )

        handler.add_method_parameter(
            "get_rho_crit",
            "set_rho_crit",
            "rho_crit",
            "density above which sink particles are created",
            default_value=(1e-16 | units.g * units.cm**-3)
        )

        handler.add_method_parameter(
            "get_r_crit",
            "set_r_crit",
            "r_crit",
            "critical radius for point mass creation"
            " (no new sinks < r_crit from existing sink)",
            default_value=(0.001 | nbody_system.length)
        )

        handler.add_method_parameter(
            "get_h_acc",
            "set_h_acc",
            "h_acc",
            "accretion radius for new sink particles",
            default_value=(0.001 | nbody_system.length)
        )

        handler.add_method_parameter(
            "get_h_soft_sinkgas",
            "set_h_soft_sinkgas",
            "h_soft_sinkgas",
            "softening length for new sink particles",
            default_value=(0.001 | nbody_system.length)
        )

        handler.add_method_parameter(
            "get_h_soft_sinksink",
            "set_h_soft_sinksink",
            "h_soft_sinksink",
            "softening length between sink particles",
            default_value=(0.001 | nbody_system.length)
        )

        handler.add_method_parameter(
            "get_f_acc",
            "set_f_acc",
            "f_acc",
            "particles < f_acc*h_acc accreted without checks",
            default_value=0.8
        )

        handler.add_method_parameter(
            "get_iexternalforce",
            "set_iexternalforce",
            "iexternalforce",
            "1=star,2=coro,3=bina,4=prdr,5=toru,6=toys,7=exte,"
            "8=spir,9=Lens,10=neut,11=Eins",
            default_value=0
        )

        handler.add_method_parameter(
            "get_irealvisc",
            "set_irealvisc",
            "irealvisc",
            "physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)",
            default_value=0
        )

        handler.add_method_parameter(
            "get_shearparam",
            "set_shearparam",
            "shearparam",
            "magnitude of shear viscosity (irealvisc=1) or alpha_SS"
            " (irealvisc=2)",
            default_value=0.1
        )

        handler.add_method_parameter(
            "get_bulkvisc",
            "set_bulkvisc",
            "bulkvisc",
            "magnitude of bulk viscosity",
            default_value=0.0
        )

        self.stopping_conditions.define_parameters(handler)

    def define_particle_sets(self, handler):
        handler.define_super_set(
            'particles',
            ['dm_particles', 'gas_particles', 'sink_particles'],
            index_to_default_set=0,
        )

        handler.define_set('dm_particles', 'index_of_the_particle')
        handler.set_new('dm_particles', 'new_dm_particle')
        handler.set_delete('dm_particles', 'delete_particle')
        handler.add_getter('dm_particles', 'get_state_dm')
        handler.add_setter('dm_particles', 'set_state_dm')
        handler.add_getter('dm_particles', 'get_mass')
        handler.add_setter('dm_particles', 'set_mass')
        handler.add_getter('dm_particles', 'get_position')
        handler.add_setter('dm_particles', 'set_position')
        handler.add_getter('dm_particles', 'get_velocity')
        handler.add_setter('dm_particles', 'set_velocity')
        handler.add_getter('dm_particles', 'get_radius')
        handler.add_setter('dm_particles', 'set_radius')

        handler.define_set('gas_particles', 'index_of_the_particle')
        handler.set_new('gas_particles', 'new_sph_particle')
        handler.set_delete('gas_particles', 'delete_particle')
        handler.add_getter('gas_particles', 'get_state_sph')
        handler.add_setter('gas_particles', 'set_state_sph')
        handler.add_getter('gas_particles', 'get_mass')
        handler.add_setter('gas_particles', 'set_mass')
        handler.add_getter('gas_particles', 'get_position')
        handler.add_setter('gas_particles', 'set_position')
        handler.add_getter('gas_particles', 'get_velocity')
        handler.add_setter('gas_particles', 'set_velocity')
        handler.add_getter('gas_particles', 'get_internal_energy')
        handler.add_setter('gas_particles', 'set_internal_energy')
        handler.add_getter('gas_particles', 'get_smoothing_length')
        handler.add_setter('gas_particles', 'set_smoothing_length')
        handler.add_getter('gas_particles', 'get_density', names=('rho',))
        handler.add_getter('gas_particles', 'get_density', names=('density',))
        handler.add_getter('gas_particles', 'get_pressure')

        handler.define_set('sink_particles', 'index_of_the_particle')
        handler.set_new('sink_particles', 'new_sink_particle')
        handler.set_delete('sink_particles', 'delete_particle')
        handler.add_getter('sink_particles', 'get_state_sink')
        handler.add_setter('sink_particles', 'set_state_sink')
        handler.add_getter('sink_particles', 'get_mass')
        handler.add_setter('sink_particles', 'set_mass')
        handler.add_getter('sink_particles', 'get_position')
        handler.add_setter('sink_particles', 'set_position')
        handler.add_getter('sink_particles', 'get_velocity')
        handler.add_setter('sink_particles', 'set_velocity')
        handler.add_getter('sink_particles', 'get_radius')
        handler.add_setter('sink_particles', 'set_radius')
        handler.add_getter('sink_particles', 'get_smoothing_length')
        handler.add_setter('sink_particles', 'set_smoothing_length')

        self.stopping_conditions.define_particle_set(handler, 'particles')

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)

        handler.add_method(
            "new_dm_particle",
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
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "new_sph_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "new_sink_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_state_dm",
            (
                handler.INDEX,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_state_dm",
            (
                handler.INDEX,
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
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_state_sph",
            (
                handler.INDEX,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_state_sph",
            (
                handler.INDEX,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.specific_energy,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_state_sink",
            (
                handler.INDEX,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                nbody_system.length,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_state_sink",
            (
                handler.INDEX,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_density",
            (
                handler.INDEX,
            ),
            (
                nbody_system.density,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_pressure",
            (
                handler.INDEX,
            ),
            (
                nbody_system.pressure,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_internal_energy",
            (
                handler.INDEX,
            ),
            (
                nbody_system.specific_energy,
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "set_internal_energy",
            (
                handler.INDEX,
                nbody_system.specific_energy,
            ),
            (
                handler.ERROR_CODE,
            )
        )


        self.stopping_conditions.define_methods(handler)
