from amuse.datamodel import Particle
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from interface import Pentacle

def test1():
    g = Pentacle(redirection="none")
    print(g.parameters)
    p = Particle()
    p.mass = 0.4 | nbody_system.mass
    p.position = [0, 1, 1.3] | nbody_system.length
    p.velocity = [-1.3, 0, -1] | nbody_system.speed
    p.radius = 0.1 | nbody_system.length
    g.particles.add_particle(p)
    print(g.particles)
    g.stop()

def test2():
    g = Pentacle(redirection="none")
    p = new_plummer_model(10)
    g.particles.add_particles(p)
    print(g.particles)
    g.stop()

def test3():
    g = Pentacle(redirection="none")
    p = new_plummer_model(10)
    g.particles.add_particles(p)
    print(g.particles)
    g.particles[0].mass = 0.01 * g.particles[0].mass
    g.particles[0].position = 0.01 * g.particles[0].position
    print(g.particles)
    g.stop()

def test4():
    # g = Pentacle(redirection="none")
    g = Pentacle()
    p = new_plummer_model(10)
    g.particles.add_particles(p)
    print(g.particles)
    g.evolve_model(0.1 | nbody_system.time)
    print(g.particles)
    g.stop()

# test3()
test4()
