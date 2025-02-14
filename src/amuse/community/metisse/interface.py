"""
Interface for metisse
"""
from amuse.community import (
    CodeInterface,
    InCodeComponentImplementation,
    LegacyFunctionSpecification,
    LiteratureReferencesMixIn,
    legacy_function,
    remote_function,
)
from amuse.community.interface import se
from amuse.datamodel import Particles, ParticlesSubset
from amuse.units import units, constants


# low level interface class
class MetisseInterface(
    CodeInterface,
    se.StellarEvolutionInterface,
    LiteratureReferencesMixIn,
):
    """
    Low level interface for metisse

    Details in publication:
        .. [#] Agrawal, P. et al. 202x
    """

    use_modules = ["metisseInterface"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self,
            name_of_the_worker="metisse_worker",
            **keyword_arguments
            )

    @remote_function
    def teststar(mass_in="d", time_in="d"):
        returns (mass_out="d", result="i")


# high level interface class
class Metisse(se.StellarEvolution):
    __interface__ = MetisseInterface

    def __init__(self, **options):
        # self.stopping_conditions = StoppingConditions(self)
        # self.stopping_conditions.supernova_detection = code.StoppingCondition('supernova_detection')
        se.StellarEvolution.__init__(self, MetisseInterface(**options), **options)

# the definition of the state model of the code
    def define_state(self, handler):
        # for example:
        # handler.set_initial_state("UNINITIALIZED")
        # handler.add_transition("!UNINITIALIZED!STOPPED", "END", "cleanup_code")
        # handler.add_transition("END", "STOPPED", "stop", False)
        # handler.add_transition(
        #     "UNINITIALIZED", "INITIALIZED", "initialize_code")
        # handler.add_method("STOPPED", "stop")
        pass

# the definition of any properties
    def define_properties(self, handler):
        # handler.add_property("name_of_the_getter", public_name="name_of_the_property")
        pass

# the definition of the parameters
    def define_parameters(self, handler):
        # handler.add_method_parameter(
        #     "name_of_the_getter",
        #     "name_of_the_setter",
        #     "parameter_name",
        #     "description",
        #     default_value = <default value>
        # )
        pass

    def define_particle_sets(self, handler):
        handler.define_set("particles", "index_of_the_star")
        handler.set_new("particles", "new_particle")
        handler.set_delete("particles", "delete_star")

        handler.add_getter("particles", "mass", "get_mass", names=("mass",))
        handler.add_getter("particles", "radius", "get_radius", names=("radius",))
        handler.add_getter("particles", "luminosity", "get_luminosity", names=("luminosity",))
        handler.add_getter("particles", "age", "get_age", names=("age",))
        handler.add_getter("particles", "stellar_type", "get_stellar_type", names=("stellar_type",))
        handler.add_getter("particles", "temperature", "get_temperature", names=("temperature",))
        handler.add_getter("particles", "time_step", "get_time_step", names=("time_step",))



class MetisseParticles(Particles):

    def __init__(self, code_interface, storage=None):
        Particles.__init__(self, storage=storage)
        self._private.code_interface = code_interface
        self.add_calculated_attribute(
            "temperature",
            self.calculate_effective_temperature,
            ["luminosity", "radius"],
        )
        self.add_function_attribute(
            "evolve_one_step", self.particleset_evolve_one_step, self.evolve_one_step
        )
        self.add_function_attribute(
            "evolve_for",
            self.particleset_evolve_for,
            self.evolve_for
        )

    def calculate_effective_temperature(self, luminosity, radius):
        return (
            (luminosity / (constants.four_pi_stefan_boltzmann * radius**2)) ** 0.25
        ).in_(
            units.K
        )

    def add_particles_to_store(self, keys, attributes=[], values=[]):
        if len(keys) == 0:
            return

        all_attributes = []
        all_attributes.extend(attributes)
        all_values = []
        all_values.extend(values)

        mapping_from_attribute_to_default_value = {
            "stellar_type": 1 | units.stellar_type,
            "radius": 0 | units.RSun,
            "luminosity": 0 | units.LSun,
            "core_mass": 0 | units.MSun,
            "CO_core_mass": 0 | units.MSun,
            "core_radius": 0 | units.RSun,
            "convective_envelope_mass": 0 | units.MSun,
            "convective_envelope_radius": 0 | units.RSun,
            "epoch": 0 | units.Myr,
            "spin": 0 | units.yr**-1,
            "main_sequence_lifetime": 0 | units.Myr,
            "age": 0 | units.Myr,
        }

        given_attributes = set(attributes)

        if "initial_mass" not in given_attributes:
            index_of_mass_attibute = attributes.index("mass")
            all_attributes.append("initial_mass")
            all_values.append(values[index_of_mass_attibute] * 1.0)

        for attribute, default_value in mapping_from_attribute_to_default_value.items():
            if attribute not in given_attributes:
                all_attributes.append(attribute)
                all_values.append(default_value.as_vector_with_length(len(keys)))

        super().add_particles_to_store(keys, all_attributes, all_values)

        added_particles = ParticlesSubset(self, keys)
        self._private.code_interface._evolve_particles(added_particles, 0 | units.yr)

    def evolve_one_step(self, particles, subset):
        self._private.code_interface._evolve_particles(
            subset.as_set(), subset.age + subset.time_step
        )

    def particleset_evolve_one_step(self, particles):
        self._private.code_interface._evolve_particles(
            particles, particles.age + particles.time_step
        )

    def evolve_for(self, particles, subset, delta_time):
        self._private.code_interface._evolve_particles(subset.as_set(), subset.age + delta_time)

    def particleset_evolve_for(self, particles, delta_time):
        self._private.code_interface._evolve_particles(particles, particles.age + delta_time)

    def get_defined_attribute_names(self):
        return ["mass", "radius"]
