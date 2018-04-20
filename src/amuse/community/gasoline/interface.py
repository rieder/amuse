from amuse.units import nbody_system
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
        # GravityFieldInterface,
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

    def new_particle(self, mass, x, y, z, vx, vy, vz, radius = 0.0):
        return self.new_dm_particle(mass, x, y, z, vx, vy, vz, radius)

    @legacy_function
    def new_dm_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.OUT)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        #function.addParameter('phi', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('radius', dtype='float64', direction=function.IN, default = 0)
        function.result_type = 'int32'
        return function

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
        function.addParameter('h_smooth', dtype='float64', direction=function.IN, default = 0)
        #function.addParameter('phi', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('metals', dtype='float64', direction=function.IN, default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_star_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.OUT)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter('tform', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('radius', dtype='float64', direction=function.IN, default = 0)
        #function.addParameter('phi', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('metals', dtype='float64', direction=function.IN, default = 0)
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

#    @legacy_function
#    def set_timestep():
#        function = LegacyFunctionSpecification()
#        function.can_handle_array = False
#        function.addParameter(
#                'timestep',
#                dtype='float64',
#                direction=function.IN)
#        #dDelta
#        function.result_type = 'int32'
#        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'radius']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state_sph():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u', 'h_smooth', 'metals']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state_star():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'tform', 'radius', 'metals']:
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
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter('radius', dtype='float64', direction=function.IN, default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state_sph():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter('h_smooth', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('metals', dtype='float64', direction=function.IN, default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state_star():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
                'index_of_the_particle',
                dtype='int32',
                direction=function.IN)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.addParameter('tform', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('radius', dtype='float64', direction=function.IN, default = 0)
        function.addParameter('metals', dtype='float64', direction=function.IN, default = 0)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['x', 'y', 'z']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['x', 'y', 'z']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_internal_energy():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_metallicity():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('metals', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def get_star_tform():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('tform', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_velocity():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['vx', 'vy', 'vz']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function    
    def set_internal_energy():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('u', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def set_metallicity():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('metals', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function    
    def set_star_tform():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('tform', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter(
                'number_of_particles',
                dtype='int32',
                direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve_model():
        function = LegacyFunctionSpecification()
        function.addParameter('time_end', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_theta():
        function = LegacyFunctionSpecification()
        function.addParameter('theta', dtype='float64', direction=function.OUT,
                description = "Barnes opening criterion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_theta():
        function = LegacyFunctionSpecification()
        function.addParameter('theta', dtype='float64', direction=function.IN,
                description = "Barnes opening criterion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_theta_2():
        function = LegacyFunctionSpecification()
        function.addParameter('theta_2', dtype='float64', direction=function.OUT,
                description = "Barnes opening criterion after a < daSwitchTheta>")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_theta_2():
        function = LegacyFunctionSpecification()
        function.addParameter('theta_2', dtype='float64', direction=function.IN,
                description = "Barnes opening criterion after a < daSwitchTheta>")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.OUT,
                description = 'fixed time step')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_time_step():
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.IN,
                description = 'fixed time step')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64', direction=function.OUT,
                description = 'time step criterion')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta():
        function = LegacyFunctionSpecification()
        function.addParameter('eta', dtype='float64', direction=function.IN,
                description = 'time step criterion')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eta_courant():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_courant', dtype='float64', direction=function.OUT,
                description = 'Courant criterion')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_eta_courant():
        function = LegacyFunctionSpecification()
        function.addParameter('eta_courant', dtype='float64', direction=function.IN,
                description = 'Courant criterion')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_alpha():
        function = LegacyFunctionSpecification()
        function.addParameter('alpha', dtype='float64', direction=function.OUT,
                description = 'Alpha constant in viscosity')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_alpha():
        function = LegacyFunctionSpecification()
        function.addParameter('alpha', dtype='float64', direction=function.IN,
                description = 'Alpha constant in viscosity')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_beta():
        function = LegacyFunctionSpecification()
        function.addParameter('beta', dtype='float64', direction=function.OUT,
                description = 'Beta constant in viscosity')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_beta():
        function = LegacyFunctionSpecification()
        function.addParameter('beta', dtype='float64', direction=function.IN,
                description = 'Beta constant in viscosity')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_nsmooth():
        function = LegacyFunctionSpecification()
        function.addParameter('nsmooth', dtype='int32', direction=function.OUT,
                description = 'number of particles to smooth over')
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_nsmooth():
        function = LegacyFunctionSpecification()
        function.addParameter('nsmooth', dtype='int32', direction=function.IN,
                description = 'number of particles to smooth over')
        function.result_type = 'int32'
        return function

# From param.dat
# gravitational collapse of an adiabatic sphere of gas.
#
# dTheta		= 0.55  # Barnes opening criterion
# dTheta2		= 0.7  # Barnes opening criterion after a < daSwitchTheta>
# dConstAlpha     = 1  # <Alpha constant in viscosity> = 1.0 or 0.5 (bBulkViscosity)
# dConstBeta      = 2  # <Beta constant in viscosity> = 2.0 or 0.5 (bBulkViscosity)
# nSmooth         = 32  # number of particles to smooth over
# iViscosityLimiter = 1  # Viscosity Limiter Type
# dEta            = 0.3  # <time step criterion> = 0.1
# dEtaCourant     = 0.2  # <Courant criterion> = 0.4

# bPeriodic	= 0  #SR Periodic boundaries
# bParaRead	= 0  #SR Ignore
# bParaWrite	= 0  #SR Ignore
# nReplicas       = 0  #SR for periodic boundaries
# achInFile	= adiabtophat.bin  #SR Ignore
# achOutName	= out  #SR Ignore
# dHubble0	= 1.0  # <dHubble0> = 0.0
# dOmega0		= 1.0  # <dOmega0> = 1.0
# dRedTo		= 0.0  # specifies final redshift for the simulation
# iCheckInterval  = 10000  # <number of timesteps between checkpoints> = 10 #SR Ignore
# dExtraStore	= 0.5  #SR ?
# bKDK		= 1  # enable/disable use of kick-drift-kick integration = +kdk
# iMaxRung	= 10  # <maximum timestep rung>
# bStandard	= 1  # output in standard TIPSY binary format = -std #SR Ignore
# bDoGravity      = 1  # enable/disable gravity (interparticle and external potentials) = +g
# bDoGas          = 1  # calculate gas/don't calculate gas = +gas
# bFastGas        = 1  # <Fast Gas Method> = 1
# dFracFastGas    = 0.3  # <Fraction of Active Particles for Fast Gas> = 0.01
# # dFracNoDomainDecomp = 0.001
# bComove         = 0  # enable/disable comoving coordinates = -cm
# bCannonical     = 1  # enable/disable use of cannonical momentum = +can
# bGeometric      = 0  # geometric/arithmetic mean to calc Grad(P/rho) = +geo
# # iGas
# bGasAdiabatic	= 1  # <Gas is Adiabatic> = +GasAdiabatic
# dConstGamma     = 1.66667  # <Ratio of specific heats> = 5/3
# dGasConst       = 0.66667  # <Gas Constant>
# dMeanMolWeight  = 1.0  # <Mean molecular weight in amu> = 1.0
# nSteps		= 3  # <number of timesteps> = 0
# dDelta          = 0.0258  # <time step>
# iOutInterval    = 1  # <number of timesteps between snapshots> = 0 #SR Ignore
# iLogInterval    = 1  # <number of timesteps between logfile outputs> = 10
# bOverwrite      = 1  # enable/disable checkpoint overwrite = -overwrite #SR Ignore
# bVDetails       = 0  # enable/disable verbose details = +vdetails
# bBulkViscosity  = 0  # <Bulk Viscosity> = 0
# bDoDensity      = 0  # enable/disable density outputs = +den
# dhMinOverSoft   = 0.25  # <Minimum h as a fraction of Softening> = 0.0


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

    def define_parameters(self, object):
        GravitationalDynamics.define_parameters(self, object)

        object.add_method_parameter(
                "get_time_step",
                "set_time_step",
                "timestep",
                "time step",
                default_value = 0.0258 | nbody_system.time,
                )

        object.add_method_parameter(
                "get_alpha",
                "set_alpha",
                "alpha",
                "Alpha constant in viscosity",
                default_value = 1.
                )

        object.add_method_parameter(
                "get_beta",
                "set_beta",
                "beta",
                "Beta constant in viscosity",
                default_value = 2.
                )
        
        object.add_method_parameter(
                "get_nsmooth",
                "set_nsmooth",
                "nsmooth",
                "number of particles to smooth over",
                default_value = 32
                )

        object.add_method_parameter(
                "get_eta",
                "set_eta",
                "eta",
                "time step criterion",
                default_value = 0.1
                )

        object.add_method_parameter(
                "get_eta_courant",
                "set_eta_courant",
                "eta_courant",
                "courant criterion",
                default_value = 0.4
                )

        object.add_method_parameter(
                "get_theta",
                "set_theta",
                "theta",
                "Barnes opening criterion",
                default_value = 0.8
                )

        object.add_method_parameter(
                "get_theta_2",
                "set_theta_2",
                "theta_2",
                "Barnes opening criterion after a < daSwitchTheta>",
                default_value = 0.55
                )

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
                "set_time_step",
                (nbody_system.time,),
                (object.ERROR_CODE,),
                )
        object.add_method(
                "set_velocity",
                (
                    object.INDEX,
                    nbody_system.speed,
                    nbody_system.speed,
                    nbody_system.speed,
                    ),
                (object.ERROR_CODE),
                )
        object.add_method(
                "get_velocity",
                (
                    object.INDEX,
                    ),
                (
                    nbody_system.speed,
                    nbody_system.speed,
                    nbody_system.speed,
                    object.ERROR_CODE
                    ),
                )
        object.add_method(
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
                    object.INDEX,
                    object.ERROR_CODE,
                    ),
                )
        object.add_method(
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
                    object.NO_UNIT,
                    ),
                (
                    object.INDEX,
                    object.ERROR_CODE,
                    ),
                )
        object.add_method(
                "new_star_particle",
                (
                    nbody_system.mass,
                    nbody_system.length,
                    nbody_system.length,
                    nbody_system.length,
                    nbody_system.speed,
                    nbody_system.speed,
                    nbody_system.speed,
                    nbody_system.time,
                    nbody_system.length,
                    object.NO_UNIT,
                    ),
                (
                    object.INDEX,
                    object.ERROR_CODE,
                    ),
                )

    def define_particle_sets(self, object):
        object.define_super_set(
                'particles',
                ['dm_particles', 'gas_particles', 'star_particles'],
                index_to_default_set=0,
                )
        
        object.define_set('dm_particles', 'index_of_the_particle')
        object.set_new('dm_particles', 'new_dm_particle')
        object.set_delete('dm_particles', 'delete_particle')
        object.add_setter('dm_particles', 'set_state')
        object.add_getter('dm_particles', 'get_state')
        object.add_setter('dm_particles', 'set_mass')
        object.add_getter('dm_particles', 'get_mass', names=('mass',))
        object.add_setter('dm_particles', 'set_position')
        object.add_getter('dm_particles', 'get_position')
        object.add_setter('dm_particles', 'set_velocity')
        object.add_getter('dm_particles', 'get_velocity')
        #object.add_getter('dm_particles', 'get_acceleration')
        object.add_setter('dm_particles', 'set_radius')
        object.add_getter('dm_particles', 'get_radius',
                names=('epsilon',))
        
        object.define_set('gas_particles', 'index_of_the_particle')
        object.set_new('gas_particles', 'new_sph_particle')
        object.set_delete('gas_particles', 'delete_particle')
        object.add_setter('gas_particles', 'set_state_sph')
        object.add_getter('gas_particles', 'get_state_sph')
        object.add_setter('gas_particles', 'set_mass')
        object.add_getter('gas_particles', 'get_mass', names = ('mass',))
        object.add_setter('gas_particles', 'set_position')
        object.add_getter('gas_particles', 'get_position')
        object.add_setter('gas_particles', 'set_velocity')
        object.add_getter('gas_particles', 'get_velocity')
        #object.add_getter('gas_particles', 'get_acceleration')
        object.add_setter('gas_particles', 'set_internal_energy')
        object.add_getter('gas_particles', 'get_internal_energy')
        object.add_setter('gas_particles', 'set_metallicity')
        object.add_getter('gas_particles', 'get_metallicity')
        object.add_getter('gas_particles', 'get_radius')
        object.add_setter('gas_particles', 'set_radius')
        #object.add_getter('gas_particles', 'get_density', names = ('rho',))
        #object.add_getter('gas_particles', 'get_density', names = ('density',))
        #object.add_getter('gas_particles', 'get_pressure')

        object.define_set('star_particles', 'index_of_the_particle')
        object.set_new('star_particles', 'new_star_particle')
        object.set_delete('star_particles', 'delete_particle')
        object.add_setter('star_particles', 'set_state_star')
        object.add_getter('star_particles', 'get_state_star')
        object.add_setter('star_particles', 'set_mass')
        object.add_getter('star_particles', 'get_mass', names = ('mass',))
        object.add_setter('star_particles', 'set_position')
        object.add_getter('star_particles', 'get_position')
        object.add_setter('star_particles', 'set_radius')
        object.add_getter('star_particles', 'get_radius')
        object.add_setter('star_particles', 'set_velocity')
        object.add_getter('star_particles', 'get_velocity')
        object.add_setter('star_particles', 'set_metallicity')
        object.add_getter('star_particles', 'get_metallicity')
        object.add_setter('star_particles', 'set_star_tform')
        object.add_getter('star_particles', 'get_star_tform')

        #self.stopping_conditions.define_particle_set(object)

