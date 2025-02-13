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
from amuse.datamodel import Particles


# low level interface class
class MetisseInterface(
    CodeInterface,
    LiteratureReferencesMixIn
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
class Metisse(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(
            self,
            MetisseInterface(**options),
            **options
        )

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

# the definition of the code data stores, either particle sets:
    def define_particle_sets(self, handler):
        # handler.define_set("particles", "index_of_the_particle")
        # handler.set_new("particles", "new_particle")
        # handler.set_delete("particles", "delete_particle")
        # handler.add_setter("particles", "set_state")
        # handler.add_getter("particles", "get_state")
        # handler.add_setter("particles", "set_mass")
        # handler.add_getter("particles", "get_mass", names=("mass",))
        pass

# and/or grids:
    def define_grids(self, handler):
        # handler.define_grid("grid",axes_names = ["x", "y"], grid_class=StructuredGrid)
        # handler.set_grid_range("grid", "_grid_range")
        # handler.add_getter("grid", "get_grid_position", names=["x", "y"])
        # handler.add_getter("grid", "get_rho", names=["density"])
        # handler.add_setter("grid", "set_rho", names=["density"])
        pass


class MetisseParticles(Particles):
    def __init__(self, code_interface, storage=None):
        Particles.__init__(self, storage=storage)
        self._private.code_interface = code_interface
