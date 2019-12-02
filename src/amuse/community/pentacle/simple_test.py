import sys
import numpy
import math
from amuse.datamodel import Particle, Particles
from amuse.units import nbody_system, units, constants
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

def test11():
    from amuse.ext.masc import new_star_cluster

    numpy.random.seed(8)
    p = new_star_cluster(
        stellar_mass=1000 | units.MSun,
        effective_radius=3.0 | units.parsec,
    )
    p.radius = 0.1 | units.parsec  # radius = switchover radius PP/PT

    M = p.mass.sum()
    R = p.position.lengths().mean()
    converter = nbody_system.nbody_to_si(M, R)
    dt = 0.01 | units.Myr
    g = Pentacle(converter)
    g.parameters.time_step = dt/1.
    print(g.parameters.time_step)
    p_in_code = g.particles.add_particles(p)
    time_end = 0.1 | units.Myr
    time = 0 | units.Myr
    channel = p_in_code.new_channel_to(p)
    step = 0
    while time < time_end:
        time += dt
        # print("Step %03i - evolving to %s" % (step, time))
        g.evolve_model(time)
        channel.copy()
        print(p[11].x.in_(units.parsec), p[11].y.in_(units.parsec))
        if (time > 9.5*dt and time < 10.5*dt):
            p_in_code.remove_particle(p_in_code[10])
        step += 1

def test12():
    from amuse.ext.masc import new_star_cluster

    numpy.random.seed(8)
    p = new_star_cluster(
        stellar_mass=1000 | units.MSun,
        effective_radius=3.0 | units.parsec,
    )
    p.radius = 0.1 | units.parsec  # radius = switchover radius PP/PT

    M = p.mass.sum()
    R = p.position.lengths().mean()
    converter = nbody_system.nbody_to_si(M, R)
    dt = 0.01 | units.Myr
    g = Pentacle(converter)
    g.parameters.time_step = dt/1.
    print(g.parameters.time_step.in_(units.Myr))
    p_in_code = g.particles.add_particles(p)
    time_end = 0.1 | units.Myr
    time = 0 | units.Myr
    channel = p_in_code.new_channel_to(p)
    step = 0
    while time < time_end:
        time += dt
        # print("Step %03i - evolving to %s" % (step, time))
        g.evolve_model(time)
        print(g.model_time.in_(units.Myr), g.particles[11].x.in_(units.parsec))
    g.stop()
    g = Pentacle(converter)
    g.parameters.time_step = dt/8.
    print(g.parameters.time_step.in_(units.Myr))
    p_in_code = g.particles.add_particles(p)
    time_end = 0.1 | units.Myr
    time = 0 | units.Myr
    channel = p_in_code.new_channel_to(p)
    step = 0
    while time < time_end:
        time += dt
        # print("Step %03i - evolving to %s" % (step, time))
        g.evolve_model(time)
        print(g.model_time.in_(units.Myr), g.particles[11].x.in_(units.parsec))
    g.stop()

def test13():
    from amuse.ext.masc import new_star_cluster

    numpy.random.seed(8)
    p = new_star_cluster(
        stellar_mass=1000 | units.MSun,
        effective_radius=3.0 | units.parsec,
    )
    p.radius = 0.1 | units.parsec  # radius = switchover radius PP/PT

    M = p.mass.sum()
    R = p.position.lengths().mean()
    converter = nbody_system.nbody_to_si(M, R)
    dt = 0.01 | units.Myr
    g = Pentacle(converter, redirection="none")
    g2 = Pentacle(converter, redirection="none")
    g.parameters.time_step = dt/1.
    g2.parameters.time_step = dt/8.
    print(g.parameters.time_step.in_(units.Myr))
    print(g2.parameters.time_step.in_(units.Myr))
    p_in_code = g.particles.add_particles(p)
    p2_in_code = g2.particles.add_particles(p)
    time_end = 0.02 | units.Myr
    time = 0 | units.Myr
    channel = p_in_code.new_channel_to(p)
    step = 0
    while time < time_end:
        time += dt
        # print("Step %03i - evolving to %s" % (step, time))
        g.evolve_model(time)
        print("g1 done 11111111111111111111111111111111111111")
        g2.evolve_model(time)
        print("g2 done 22222222222222222222222222222222222222")
        # print("g dt:   ", g.model_time.in_(units.Myr), g.particles[11].x.in_(units.parsec))
        # print("g dt/8: ", g2.model_time.in_(units.Myr), g2.particles[11].x.in_(units.parsec))
        print(
            (
                g.particles.position
                - g2.particles.position
            ).lengths().mean().in_(units.parsec)
        )
    g.stop()
    g2.stop()

def test14():
    numpy.random.seed(8)
    from amuse.lab import new_solar_system
    p = new_solar_system()
    p.radius = 9 | units.AU  # radius = switchover radius PP/PT

    M = p.mass.sum()
    R = p.position.lengths().mean()
    converter = nbody_system.nbody_to_si(M, R)
    dt = 0.01 | units.yr
    g = Pentacle(converter)
    g.parameters.time_step = dt/1.
    print(g.parameters.time_step.in_(units.Myr))
    p_in_code = g.particles.add_particles(p)
    time_end = 2 | units.yr
    time = 0 | units.yr
    channel = p_in_code.new_channel_to(p)
    step = 0
    while time < time_end:
        time += dt
        # print("Step %03i - evolving to %s" % (step, time))
        g.evolve_model(time)
        sun = g.particles[0]
        print(
            "%s %s %s %s %s" % (
                g.model_time,
                (
                    g.particles[1].position - sun.position
                ).lengths().in_(units.au),
                (
                    g.particles[3].position - sun.position
                ).lengths().in_(units.au),
                (
                    g.particles[5].position - sun.position
                ).lengths().in_(units.au),
                (
                    g.particles[7].position - sun.position
                ).lengths().in_(units.au),
            )
        )
        # print("g dt:   ", g.model_time.in_(units.Myr), g.particles[11].x.in_(units.parsec))
    g.stop()

def test15():
    # It seems there is an error with (very) small time_step values?
    # In this example, dt/4 will update particles at every step while dt/8 will
    # not update ~half of the particles each timestep...
    # Apparently dt_soft should not be smaller than ~0.00037 ??
    from amuse.ext.masc import new_star_cluster

    numpy.random.seed(8)
    p = new_star_cluster(
        stellar_mass=1000 | units.MSun,
        effective_radius=5.0 | units.parsec,
    )
    p.radius = 0.1 | units.parsec  # radius = switchover radius PP/PT
    p = new_plummer_model(65536)
    p.radius = (2/256) | nbody_system.length
    M = p.mass.sum()
    R = p.position.lengths().mean()
    # converter = nbody_system.nbody_to_si(M, R)
    # print(converter.to_nbody(0.1 | units.parsec))
    # exit()
    dt = 0.0625 | nbody_system.time  # | units.Myr
    g = Pentacle(number_of_workers=1, redirection="none")
    # print(1/len(p))
    # g.parameters.time_step = 0.00037 | nbody_system.time
    g.parameters.time_step = 1./256. | nbody_system.time
    # p.radius = converter.to_si(
    #     g.parameters.time_step * (2 | nbody_system.speed)
    # )
    # print(p[0].radius.in_(units.parsec))
    # exit()
    # g.parameters.time_step = dt/10.
    # print(converter.to_nbody(g.parameters.time_step))
    # print(g.parameters.time_step.in_(units.Myr))
    # exit()
    p_in_code = g.particles.add_particles(p)
    time_end = 16*g.parameters.time_step  # 0.05 | units.Myr
    time = 0 | nbody_system.time
    channel = p_in_code.new_channel_to(p)
    step = 0
    initial_pos = g.particles.x
    while time < time_end:
        time += g.parameters.time_step
        print("Step %03i - evolving to %s" % (step, time))
        previous_pos = g.particles.x
        g.evolve_model(time)
        current_pos = g.particles.x
        print(
            len(
                p[
                    # current_pos == initial_pos
                    current_pos == previous_pos
                ]
            )
        )
        # print(g.model_time.in_(units.Myr), g.particles[11].x.in_(units.parsec))
        step += 1
    g.stop()

def test16():
    # It seems there is an error with (very) small time_step values?
    # In this example, dt/4 will update particles at every step while dt/8 will
    # not update ~half of the particles each timestep...
    # Apparently dt_soft should not be smaller than ~0.00037 ??
    from amuse.ext.masc import new_star_cluster

    numpy.random.seed(8)
    p = new_plummer_model(2)

    dt = (1.0/4096) | nbody_system.time
    g = Pentacle(number_of_workers=1, redirection="none")
    # print(1/len(p))
    # g.parameters.time_step = 0.00037 | nbody_system.time
    g.parameters.time_step = dt
    p.radius = g.parameters.time_step * (2 | nbody_system.speed)
    p_in_code = g.particles.add_particles(p)
    print("Added particles")
    time_end = 50*g.parameters.time_step
    time = 0 | nbody_system.time
    channel = p_in_code.new_channel_to(p)
    step = 0
    initial_pos = g.particles.x
    while time < time_end:
        time += g.parameters.time_step
        previous_pos = g.particles.x
        g.evolve_model(time)
        current_pos = g.particles.x
        print(
            len(
                p[
                    # current_pos == initial_pos
                    current_pos == previous_pos
                ]
            )
        )
        # print(g.model_time.in_(units.Myr), g.particles[11].x.in_(units.parsec))
    g.stop()

def test17():
    particles = Particles(2)
    particles.mass = 1.0 | units.MSun
    particles.radius = 10.0 | units.AU
    particles.position = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] | units.AU
    particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
    particles[1].vy = (constants.G * (2.0 | units.MSun) / (2.0 | units.AU)).sqrt()
    particles.move_to_center()
    
    converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
    instance = Pentacle(converter, redirection="none")
    instance.parameters.time_step = 8 * 0.0125 * math.pi * particles[0].x / particles[0].vy
    instance.commit_parameters()
    instance.particles.add_particles(particles)
    instance.commit_particles()
    primary = instance.particles[0]
    
    P = 2 * math.pi * primary.x / primary.vy
    # print(P / instance.parameters.time_step)
    
    position_at_start = primary.position.x
    dt = instance.parameters.time_step
    time = 0 * dt
    while time < (dt):
        time += instance.parameters.time_step
        instance.evolve_model(time)
        print(time/P)
    # instance.evolve_model(P / 4.0)
    # print("at P/4")
    
    # instance.evolve_model(P / 2.0)
    
    # instance.evolve_model(P)
    
    instance.stop()


# testDefaultParameters()
# test3()
test17()
