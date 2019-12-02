import math

from amuse.test.amusetest import TestWithMPI

from amuse.datamodel import Particles
from amuse.units import units, constants, nbody_system
from amuse.ic.plummer import new_plummer_model

from .interface import PentacleInterface
from .interface import Pentacle


default_options = dict() # redirection="none")

class PentacleInterfaceTests(TestWithMPI):
    
    def test_new_instance(self):
        instance = self.new_instance_of_an_optional_code(
            PentacleInterface,
            **default_options
        )
        self.assertEqual(0, instance.initialize_code())
        instance.stop()

    def test_new_particles(self):
        instance = self.new_instance_of_an_optional_code(
            PentacleInterface,
            **default_options
        )
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())

        index, error = instance.new_particle(
            mass=11.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0, radius=0.1
        )
        self.assertEqual(0, error)
        self.assertEqual(1, index)
        index, error = instance.new_particle(
            mass=21.0, x=10.0, y=0.0, z=0.0, vx=10.0, vy=0.0, vz=0.0, radius=0.2
        )
        self.assertEqual(0, error)
        self.assertEqual(2, index)
        self.assertEqual(0, instance.commit_particles())
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
        self.assertEqual(0, retrieved_state1['__result'])
        self.assertEqual(0, retrieved_state2['__result'])
        self.assertEqual(11.0, retrieved_state1['mass'])
        self.assertEqual(21.0, retrieved_state2['mass'])
        self.assertEqual(0.0, retrieved_state1['x'])
        self.assertEqual(10.0, retrieved_state2['x'])

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test_particle_property_setters_getters(self):
        instance = self.new_instance_of_an_optional_code(
            PentacleInterface,
            **default_options
        )
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual([1, 0], list(instance.new_particle(0.01,  1, 0, 0,  0, 1, 0, 0.1).values()))
        self.assertEqual([2, 0], list(instance.new_particle(0.02, -1, 0, 0,  0,-1, 0, 0.1).values()))
        self.assertEqual(0, instance.commit_particles())
        
        # getters
        mass, result = instance.get_mass(1)
        self.assertAlmostEqual(0.01, mass)
        self.assertEqual(0,result)
        radius, result = instance.get_radius(2)
        self.assertAlmostEqual(0.1, radius)
        self.assertEqual(0,result)
        self.assertEqual(-3, instance.get_mass(3)['__result']) # Particle not found
        self.assertEqual([ 1, 0, 0,  0], list(instance.get_position(1).values()))
        self.assertEqual([-1, 0, 0,  0], list(instance.get_position(2).values()))
        self.assertEqual([ 0, 1, 0,  0], list(instance.get_velocity(1).values()))
        self.assertEqual([ 0,-1, 0,  0], list(instance.get_velocity(2).values()))
        
        # setters
        self.assertEqual(0, instance.set_state(1, 0.01, 1,2,3, 4,5,6, 0.1))
        self.assertEqual([0.01, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], list(instance.get_state(1).values()))
        self.assertEqual(0, instance.set_mass(1, 0.02))
        self.assertEqual([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], list(instance.get_state(1).values()))
        self.assertEqual(0, instance.set_radius(1, 0.2))
        self.assertEqual([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.2, 0], list(instance.get_state(1).values()))
        self.assertEqual(0, instance.set_position(1, 10,20,30))
        self.assertEqual([0.02, 10.0,20.0,30.0, 4.0,5.0,6.0, 0.2, 0], list(instance.get_state(1).values()))
        self.assertEqual(0, instance.set_velocity(1, 40,50,60))
        self.assertEqual([0.02, 10.0,20.0,30.0, 40.0,50.0,60.0, 0.2, 0], list(instance.get_state(1).values()))

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test_parameters(self):
        instance = self.new_instance_of_an_optional_code(PentacleInterface, **default_options)
        self.assertEqual(0, instance.initialize_code())
        
        self.assertEqual([1.0e-2, 0], list(instance.get_eps2().values()))
        
        self.assertEqual(0, instance.set_eps2(0.2))
        self.assertEqual([0.2, 0], list(instance.get_eps2().values()))
        
        self.assertEqual([0.1, 0], list(instance.get_eta().values()))
        
        self.assertEqual(0, instance.set_eta(0.01))
        self.assertEqual([0.01, 0], list(instance.get_eta().values()))
        
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

class PentacleTests(TestWithMPI):

    default_converter = nbody_system.nbody_to_si(1.0e4 | units.MSun, 1.0 | units.AU)

    def new_sun_earth_system(self):
        particles = Particles(2)
        particles.mass = [1.0, 3.0037e-6] | units.MSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * particles.total_mass() / (1.0 | units.AU)).sqrt()
        return particles

    def test_initialization(self):
        instance = self.new_instance_of_an_optional_code(Pentacle, self.default_converter, **default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()

    def test_parameters(self):
        instance = self.new_instance_of_an_optional_code(Pentacle, self.default_converter, **default_options)
        instance.initialize_code()
        
        self.assertEqual(
            instance.parameters.epsilon_squared, 
            instance.unit_converter.to_si(1.0e-2 | nbody_system.length**2)
        )
        self.assertEqual(instance.parameters.time_step_parameter, 0.1)
        self.assertEqual(
            instance.parameters.time_step,
            instance.unit_converter.to_si(1/256 | nbody_system.time)
        )
        
        # for par, value in [
        #         # ('opening_angle', 0.4),
        # ]:
        #     self.assertEqual(instance.unit_converter.to_si(value), 
        #         getattr(instance.parameters, par))
        #         
        #     if hasattr(value, 'unit'):
        #         new_value = 3.0 | value.unit
        #     else:
        #         new_value = 3.0
        #         
        #     setattr(instance.parameters, par, new_value)
        #     self.assertEqual(instance.unit_converter.to_si(new_value),
        #         getattr(instance.parameters, par))
        
        instance.commit_parameters()
        p = instance.parameters
        
        instance.stop()

    def test_instance(self):
        instance = self.new_instance_of_an_optional_code(Pentacle, self.default_converter, **default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()
        
        self.assertAlmostEqual(instance.particles.mass, [1.0, 3.0037e-6] | units.MSun)
        self.assertAlmostEqual(instance.particles.position, 
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU)
        self.assertAlmostEqual(instance.particles.velocity, 
            [[0.0, 0.0, 0.0], [0.0, 29.7885, 0.0]] | units.km / units.s, 3)
        
        instance.cleanup_code()
        instance.stop()

    def test_evolve_oneparticle(self):
        particles = Particles(1)
        particles.mass = 1.0 | units.MSun
        # instance.parameters.rcut_out_star_star = 10.0 | units.AU
        # particles.radius = 10.0 | units.AU
        particles.position = [[0.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[1.0, 0.0, 0.0]] | units.AU / units.yr
        
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Pentacle(converter) #self.new_instance_of_an_optional_code(Pentacle, converter, **default_options)
        instance.initialize_code()
        instance.parameters.time_step = 0.01 | units.yr
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]
        
        P = 0.1 | units.yr
        
        position_at_start = primary.position.x
        dt = instance.parameters.time_step
        time = 0 * dt
        
        instance.evolve_model(P / 2.0)
        self.assertAlmostRelativeEqual(0.05 | units.AU, primary.position.x, 3)
        
        # instance.evolve_model(P)
        # self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)
        
        instance.cleanup_code()
        instance.stop()

    def test_evolve_twoparticles_a(self):
        particles = Particles(2)
        particles.mass = [1.0, 0.001] | units.MSun
        # instance.parameters.rcut_out_star_star = 10.0 | units.AU
        # particles.radius = 10.0 | units.AU
        particles.position = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.AU / units.yr
        
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Pentacle(converter) #self.new_instance_of_an_optional_code(Pentacle, converter, **default_options)
        instance.initialize_code()
        instance.parameters.time_step = 0.01 | units.yr
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]
        
        P = 0.1 | units.yr
        
        position_at_start = primary.position.x
        dt = instance.parameters.time_step
        time = 0 * dt
        
        instance.evolve_model(P)
        self.assertAlmostRelativeEqual(primary.position.x, 0.1 | units.AU, 3)
        
        # instance.evolve_model(P)
        # self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)
        
        instance.cleanup_code()
        instance.stop()

    def test_evolve_twoparticles(self):
        # from amuse.community.ph4.interface import ph4
        # from amuse.community.bhtree.interface import BHTree
        particles = Particles(2)
        particles.mass = 1.0 | units.MSun
        # instance.parameters.rcut_out_star_star = 10.0 | units.AU
        particles.radius = 0.1 | units.AU
        particles.position = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * (2.0 | units.MSun) / (2.0 | units.AU)).sqrt()
        particles.move_to_center()
        
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Pentacle(converter) #self.new_instance_of_an_optional_code(Pentacle, converter, **default_options)
        instance.initialize_code()
        instance.parameters.time_step = 0.0125 * math.pi * particles[0].x / particles[0].vy
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]
        
        P = 2 * math.pi * primary.x / primary.vy
        
        position_at_start = primary.position.x
        dt = instance.parameters.time_step
        time = 0 * dt
        instance.evolve_model(P / 4.0)
        #self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 2)
        
        instance.evolve_model(P / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 2)
        
        instance.evolve_model(P)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)
        
        instance.cleanup_code()
        instance.stop()

    def test_earthsun(self):
        # from amuse.community.ph4.interface import ph4
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = self.new_instance_of_an_optional_code(Pentacle, converter, **default_options)
        instance.initialize_code()
        instance.parameters.time_step = 0.005 | units.yr
        instance.commit_parameters()
        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()
        earth = instance.particles[1]
        
        position_at_start = earth.position.x
        instance.evolve_model(0.25 | units.yr)
        #self.assertAlmostRelativeEqual(position_at_start, earth.position.y, 2)
        
        instance.evolve_model(0.5 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, -earth.position.x, 2)
        
        instance.evolve_model(1.0 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, earth.position.x, 2)
        
        instance.cleanup_code()
        instance.stop()

