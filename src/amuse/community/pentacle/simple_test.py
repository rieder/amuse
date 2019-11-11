import sys
import numpy
from amuse.datamodel import Particle
from amuse.units import nbody_system, units
from amuse.ic.plummer import new_plummer_model
from amuse.io import read_set_from_file
from interface import Pentacle

def testDefaultParameters():
    g = Pentacle()
    print(g.parameters)
    g.stop()

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
    p = new_plummer_model(1000)
    g.particles.add_particles(p)
    # print(g.particles)
    for i in range(10):
        time = (i/10.) | nbody_system.time
        g.evolve_model(time)
        print(time, g.particles[0].position)
    # print(g.particles)
    g.stop()

def test5():
    numpy.random.seed(5)
    g = Pentacle(redirection="none")
    p = new_plummer_model(1024)
    # p = read_set_from_file(sys.argv[1], "amuse")
    g.parameters.epsilon_squared = ((4/len(p)) | nbody_system.length)**2
    g.particles.add_particles(p)
    dt = (1.0 / 256.) | nbody_system.time
    t_end = 0.125 | nbody_system.time
    dt = t_end
    time = 0 | nbody_system.time
    while time < t_end:
        time += dt
        g.evolve_model(time)
        print(time, g.particles[0].position)
    # print(g.particles)
    g.stop()

def test6():
    numpy.random.seed(5)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from amuse.plot import plot
    g = Pentacle(redirection="none")
    # g = Pentacle(number_of_workers=1)
    # p = new_plummer_model(1024)
    p = read_set_from_file(sys.argv[1], "amuse")
    g.parameters.epsilon_squared = ((4/len(p)) | nbody_system.length)**2
    g.particles.add_particles(p)
    # print(g.parameters)
    # exit()
    t = [] | nbody_system.time
    x = [] | nbody_system.length
    y = [] | nbody_system.length
    z = [] | nbody_system.length
    dt = (1.0 / 256.) | nbody_system.time
    t_end = 1.5 | nbody_system.time
    # dt = t_end
    time = 0 | nbody_system.time
    while time < t_end:
        time += dt
        g.evolve_model(time)
        t.append(time)
        x.append(g.particles[0].x)
        y.append(g.particles[0].y)
        z.append(g.particles[0].z)
        print(time, g.particles[0].position)
    # print(g.particles)
    g.stop()
    plot(t, x)
    plot(t, y)
    plot(t, z)
    plt.savefig("test6-xyz.png")

def test7():
    numpy.random.seed(5)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from amuse.plot import plot
    from amuse.io import write_set_to_file
    g = Pentacle(redirection="none")
    # g = Pentacle(number_of_workers=1)
    N = 8192
    p = new_plummer_model(N)

    q = new_plummer_model(N)

    p.x -= 2 | nbody_system.length
    q.x += 2 | nbody_system.length
    p.vy += 4 | nbody_system.speed
    q.vy -= 4 | nbody_system.speed
    r = Particle()
    r.position = [0,0,0] | nbody_system.length
    r.velocity = [0,0,0] | nbody_system.speed
    r.mass = 40 | nbody_system.mass
    g.parameters.epsilon_squared = ((4/len(p)) | nbody_system.length)**2
    g.particles.add_particles(p)
    g.particles.add_particles(q)
    g.particles.add_particle(r)
    # print(g.parameters)
    # exit()
    t = [] | nbody_system.time
    x = [] | nbody_system.length
    y = [] | nbody_system.length
    z = [] | nbody_system.length
    x2 = [] | nbody_system.length
    y2 = [] | nbody_system.length
    z2 = [] | nbody_system.length
    dt = (1.0 / 128.) | nbody_system.time
    t_end = 1.0 | nbody_system.time
    # dt = t_end
    time = 0 | nbody_system.time
    while time < t_end:
        time += dt
        g.evolve_model(time)
        t.append(time)
        x.append(g.particles[0].x)
        y.append(g.particles[0].y)
        z.append(g.particles[0].z)
        x2.append(g.particles[N].x)
        y2.append(g.particles[N].y)
        z2.append(g.particles[N].z)
        print(time, g.particles[0].position)
        print(time, g.particles[N].position)
    # print(g.particles)

    plot(x, y)
    plot(x2, y2)
    plt.savefig("test7-xy.png")
    write_set_to_file(g.particles, "test7-final.amuse", "amuse")
    g.stop()

def test10():
    from amuse.ext.masc import new_star_cluster
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    numpy.random.seed(8)
    p = new_star_cluster(
        stellar_mass=1000 | units.MSun,
        effective_radius=3.0 | units.parsec,
    )
    p.radius = 0.1 | units.parsec  # radius = switchover radius PP/PT
    q = new_star_cluster(
        stellar_mass=1000 | units.MSun,
        effective_radius=3.0 | units.parsec,
    )
    q.radius = 0.1 | units.parsec

    p.x -= 3.5 | units.parsec
    q.x += 3.5 | units.parsec
    p.vx += 0.5 | units.kms
    q.vy -= 1 | units.kms
    M = p.mass.sum()
    R = p.position.lengths().mean()
    converter = nbody_system.nbody_to_si(M, R)
    dt = 0.05 | units.Myr
    g = Pentacle(converter)
    g.parameters.time_step = dt/5
    p_in_code = g.particles.add_particles(p)
    time_end = 10 | units.Myr
    time = 0 | units.Myr
    channel = p_in_code.new_channel_to(p)
    step = 0
    second_added = False
    while time < time_end:
        time += dt
        if not second_added:
            if time > 0.5 | units.Myr:
                q_in_code = g.particles.add_particles(q)
                channelq = q_in_code.new_channel_to(q)
                second_added = True
        print("Step %03i - evolving to %s" % (step, time))
        g.evolve_model(time)
        channel.copy()
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, aspect=1)
        ax.scatter(
            p.x.value_in(units.parsec),
            p.y.value_in(units.parsec),
            s=numpy.sqrt(p.mass.value_in(units.MSun)),
            edgecolors="none",
        )
        if second_added:
            channelq.copy()
            ax.scatter(
                q.x.value_in(units.parsec),
                q.y.value_in(units.parsec),
                s=numpy.sqrt(q.mass.value_in(units.MSun)),
                edgecolors="none",
            )
        ax.set_xlim((-10, 10))
        ax.set_ylim((-10, 10))
        plt.savefig("test10-%03i.png"%step, dpi=150)
        plt.close(fig)
        print("COMV (cluster 1): ")
        print(p.center_of_mass_velocity().in_(units.kms))
        step += 1


# testDefaultParameters()
# test3()
test10()
