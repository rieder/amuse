"""
Interface for metisse
"""

from amuse.community import (
    CodeInterface,
    LiteratureReferencesMixIn,
    # LegacyFunctionSpecification,
    # legacy_function,
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
    Low level interface for METISSE

    Details in publication:
        .. [#] Agrawal, P. et al., 2020, https://doi.org/10.1093/mnras/staa2264
        .. [#] Agrawal, P. et al., 2023, https://doi.org/10.1093/mnras/stad2334
    """

    use_modules = ["metisseInterface"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="metisse_worker", **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)
        self.model_time = 0.0 | units.julianyr

    # Remote functions - getters and setters
    # Note that we should maybe use SI units rather than derived (MSun etc), at
    # least while these are not certain to be the same in the code and in
    # AMUSE...
    @remote_function(can_handle_array=True)
    def get_initial_mass(index_of_the_star="i"):
        returns (mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_epoch(index_of_the_star="i"):
        returns (epoch="d" | units.julianyr)

    @remote_function(can_handle_array=True)
    def get_core_mass(index_of_the_star="i"):
        returns (core_mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_core_radius(index_of_the_star="i"):
        returns (core_radius="d" | units.RSun)

    @remote_function(can_handle_array=True)
    def get_convective_envelope_mass(index_of_the_star="i"):
        returns (convective_envelope_mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_convective_envelope_radius(index_of_the_star="i"):
        returns (convective_envelope_radius="d" | units.RSun)

    @remote_function(can_handle_array=True)
    def get_CO_core_mass(index_of_the_star="i"):
        returns (CO_core_mass="d" | units.MSun)

    @remote_function(can_handle_array=True)
    def get_main_sequence_lifetime(index_of_the_star="i"):
        returns (main_sequence_lifetime="d" | units.Myr)

    @remote_function(must_handle_array=True)
    def evolve_stars(index_of_the_star="i", time_delta="d" | units.Myr):
        returns (error="i")

    # getters and setters for tracks
    # metallicity_dir (string)
    # metallicity_dir_he (string)
    # z_accuracy_limit (float)
    # mass_accuracy_limit (float)

    @remote_function
    def get_metallicity_dir():
        returns (metallicity_dir="s")

    @remote_function
    def set_metallicity_dir(metallicity_dir="s"):
        returns ()

    @remote_function
    def get_metallicity_dir_he():
        returns (metallicity_dir_he="s")

    @remote_function
    def set_metallicity_dir_he(metallicity_dir_he="s"):
        returns ()

    @remote_function
    def get_z_accuracy_limit():
        returns (z_accuracy_limit="d")

    @remote_function
    def set_z_accuracy_limit(z_accuracy_limit="d"):
        returns ()

    @remote_function
    def get_mass_accuracy_limit():
        returns (mass_accuracy_limit="d")

    @remote_function
    def set_mass_accuracy_limit(mass_accuracy_limit="d"):
        returns ()

    # getters and setters for miscellaneous controls
    # verbose (bool)
    # construct_postagb_track (bool)

    @remote_function
    def get_verbose():
        returns (verbose="b")

    @remote_function
    def set_verbose(verbose="b"):
        returns ()

    @remote_function
    def get_construct_postagb_track():
        returns (construct_postagb_track="b")

    @remote_function
    def set_construct_postagb_track(construct_postagb_track="b"):
        returns ()

    # getters and setters for parameters
    # initial_metallicity(real)
    # wd_mass_scheme (string, 256)
    # use_initial_final_mass_relation(bool)
    # bhns_mass_scheme (string, 256)
    # max_ns_mass (real)
    # allow_electron_capture (bool)

    @remote_function
    def get_initial_metallicity():
        returns (initial_metallicity="d")

    @remote_function
    def set_initial_metallicity(initial_metallicity="d"):
        returns ()

    @remote_function
    def get_wd_mass_scheme():
        returns (wd_mass_scheme="s")

    @remote_function
    def set_wd_mass_scheme(wd_mass_scheme="s"):
        returns ()

    @remote_function
    def get_use_initial_final_mass_relation():
        returns (use_initial_final_mass_relation="b")

    @remote_function
    def set_use_initial_final_mass_relation(use_initial_final_mass_relation="b"):
        returns ()

    @remote_function
    def get_bhns_mass_scheme():
        returns (bhns_mass_scheme="s")

    @remote_function
    def set_bhns_mass_scheme(bhns_mass_scheme="s"):
        returns ()

    @remote_function
    def get_max_ns_mass():
        returns (max_ns_mass="d")

    @remote_function
    def set_max_ns_mass(max_ns_mass="d"):
        returns ()

    @remote_function
    def get_allow_electron_capture():
        returns (allow_electron_capture="b")

    @remote_function
    def set_allow_electron_capture(allow_electron_capture="b"):
        returns ()


# high level interface class
class Metisse(se.StellarEvolution):
    """
    High level interface for METISSE
    """
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

        

        # Track parameters
        handler.add_method_parameter(
            "get_metallicity_dir",
            "set_metallicity_dir",
            "metallicity_dir",
            "Location of the tracks",
            default_value="./",
        )

        handler.add_method_parameter(
            "get_metallicity_dir_he",
            "set_metallicity_dir_he",
            "metallicity_dir_he",
            "Location of the He tracks",
            default_value="./",
        )

        handler.add_method_parameter(
            "get_z_accuracy_limit",
            "set_z_accuracy_limit",
            "z_accuracy_limit",
            "Metallicity accuracy limit",
            default_value=1.0e-2,
        )

        handler.add_method_parameter(
            "get_mass_accuracy_limit",
            "set_mass_accuracy_limit",
            "mass_accuracy_limit",
            "Mass accuracy limit",
            default_value=1.0e-4,
        )

        # handlers for parameters:
        # initial_metallicity
        # wd_mass_scheme
        # use_initial_final_mass_relation
        # bhns_mass_scheme
        # max_ns_mass
        # allow_electron_capture

        handler.add_method_parameter(
            "get_initial_metallicity",
            "set_initial_metallicity",
            "initial_metallicity",
            "Initial metallicity",
            default_value=-1.0,
        )

        handler.add_method_parameter(
            "get_wd_mass_scheme",
            "set_wd_mass_scheme",
            "wd_mass_scheme",
            (
                "White Dwarf (WD) luminosity calculation method:\n"
                "(1) \"Mestel\" - Shapiro S. L., Teukolsky S. A., 1983\n"
                "(2) \"Modified_mestel\" - Hurley J. R., Shara M. M., 2003"
            ),
            default_value="Modified_mestel",
        )

        handler.add_method_parameter(
            "get_use_initial_final_mass_relation",
            "set_use_initial_final_mass_relation",
            "use_initial_final_mass_relation",
            (
                "If True use the initial final mass relation for white dwarfs "
                "from Han, Z., Posialowski, P., Eggleton, P. P., 1995."
            ),
            default_value=False,
        )

        handler.add_method_parameter(
            "get_bhns_mass_scheme",
            "set_bhns_mass_scheme",
            "bhns_mass_scheme",
            (
                "Neutron Star/Black Hole (NS/BH) type and mass calculation method:\n"
                "(1) \"original_SSE\" - Hurley et al. 2000\n"
                "(2) \"Belczynski2002\" - Belczynski et al. 2002\n"
                "(3) \"Belczynski2008\" - Belczynski et al. 2008\n"
                "(4) \"Eldridge_Tout2004\" - Eldridge J. J., Tout C. A., 2004"
            ),
            default_value="Belczynski2008",
        )

    def define_particle_sets(self, handler):
        handler.define_set("particles", "index_of_the_star")
        handler.set_new("particles", "new_particle")
        handler.set_delete("particles", "delete_star")

        handler.add_getter("particles", "get_mass", names=("mass",))
        handler.add_getter("particles", "get_radius", names=("radius",))
        handler.add_getter("particles", "get_age", names=("age",))
        handler.add_getter(
            "particles", "get_time_step", names=("time_step",)
        )
        handler.add_getter(
            "particles", "get_temperature", names=("temperature",)
        )
        handler.add_getter(
            "particles", "get_luminosity", names=("luminosity",)
        )
        handler.add_getter(
            "particles", "get_stellar_type", names=("stellar_type",)
        )
        # handler.add_getter("particles", "get_spin", names=("spin",))
        handler.add_getter("particles", "get_epoch", names=("epoch",))
        handler.add_getter(
            "particles",
            "get_main_sequence_lifetime",
            names=("main_sequence_lifetime",),
        )
        handler.add_getter(
            "particles", "get_core_mass", names=("core_mass",)
        )
        handler.add_getter(
            "particles", "get_CO_core_mass", names=("CO_core_mass",)
        )
        handler.add_getter(
            "particles", "get_core_radius", names=("core_radius",)
        )
        handler.add_getter(
            "particles",
            "get_convective_envelope_mass",
            names=("convective_envelope_mass",),
        )
        handler.add_getter(
            "particles",
            "get_convective_envelope_radius",
            names=("convective_envelope_radius",),
        )
        handler.add_getter(
            "particles", "get_initial_mass", names=("initial_mass",)
        )

    def evolve_model(self, end_time=None, keep_synchronous=True):
        if not keep_synchronous:
            self._evolve_particles(self.particles, self.particles.time_step + self.particles.age)
            return

        if end_time is None:
            end_time = self.model_time + min(self.particles.time_step)
        self.evolve_stars(self.particles, end_time - self.model_time + self.particles.age)
        self.model_time = end_time


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
            "evolve_for", self.particleset_evolve_for, self.evolve_for
        )

    def calculate_effective_temperature(self, luminosity, radius):
        return (
            (luminosity / (constants.four_pi_stefan_boltzmann * radius**2)) ** 0.25
        ).in_(units.K)

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
        self._private.code_interface._evolve_particles(
            subset.as_set(), subset.age + delta_time
        )

    def particleset_evolve_for(self, particles, delta_time):
        self._private.code_interface._evolve_particles(
            particles, particles.age + delta_time
        )

    def get_defined_attribute_names(self):
        return ["mass", "radius"]
