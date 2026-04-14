import math
import numpy as np
from amuse.lab import *
from amuse.couple import bridge
from plummer_potential import Plummer_potential
#from paczynski_potential import Paczynski_potential
#from hermite.interface import Hermite

from matplotlib import pyplot as plt

# #BOOKLISTSTART1# #
class CodeWithFriction(bridge.GravityCodeInField):
    
    def kick_with_field_code(self, particles, field_code, dt):
        self.LnL = 3.7

        R = particles.position.length()
        vx = particles.vx.mean()
        vy = particles.vy.mean()
        vz = particles.vz.mean()
        rho = field_code.density(R)
        vc = field_code.circular_velocity(R)
        X = 0.34 
        m = particles.mass.sum()
        ax = -4*np.pi*self.LnL*constants.G**2 * rho*m*(vx/vc**3)*X
        ay = -4*np.pi*self.LnL*constants.G**2 * rho*m*(vy/vc**3)*X
        az = -4*np.pi*self.LnL*constants.G**2 * rho*m*(vz/vc**3)*X
        self.update_velocities(particles, dt, ax, ay, az)

    def drift(self, tend): 
        pass

def get_system_state_relative_to_SMBH(time, black_holes):
    t = [] | units.Myr
    r = [] | units.pc
    SMBH = black_holes[black_holes.name=="SMBH"]
    for bi in black_holes-SMBH:
        t.append(time)
        r.append((bi.position-SMBH.position).length())
    return t, r

def get_system_state(time, black_holes):
    t = [] | units.Myr
    r = [] | units.pc
    for bi in black_holes:
        t.append(time)
        r.append(bi.position.length())
    return t, r
# #BOOKLISTSTOP1# #

def evolve_cluster_in_potential(gravity, cluster_gravity,
                                black_holes, t_end, dt,
                                channel_to_framework):

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    x = []
    y = []
    t, r = get_system_state(gravity.model_time, black_holes)
    x.append(t)
    y.append(r)
    while gravity.model_time < t_end:

        gravity.evolve_model(gravity.model_time + dt)
        channel_to_framework.copy()

        t, r = get_system_state(gravity.model_time, black_holes)
        x.append(t)
        y.append(r)        
#        x.append(gravity.model_time.value_in(units.Myr))
#        y.append(black_holes.position.length().value_in(units.parsec))

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy
        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print("T=", gravity.model_time.in_(units.Myr), end=' ') 
        print("E= ", Etot.in_(units.erg), "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot) 
        Etot_prev = Etot
    return x, y

def integrate_black_holes_in_potential(black_holes, potential, t_end, dt, dt_bridge, converter):

    cluster_gravity = Hermite(converter)
    #cluster_gravity.parameters.light_speed = 1*constants.c

    cluster_gravity.particles.add_particles(black_holes)
    channel_from_gravity_to_framework = cluster_gravity.particles.new_channel_to(black_holes)

    friction_code = CodeWithFriction(cluster_gravity, (potential,), do_sync=True, verbose=False, radius_is_eps=False, h_smooth_is_eps=False, zero_smoothing=False)
        
    gravity = bridge.Bridge(use_threading=False)
    gravity.add_system(cluster_gravity, (potential,) )
    gravity.add_code(friction_code)
    gravity.timestep = dt_bridge

    x, y = evolve_cluster_in_potential(gravity, cluster_gravity,
                                       black_holes, t_end, dt,
                                       channel_from_gravity_to_framework)
    gravity.stop()
    return x, y

def plot_orbit(x, y):
    x = np.array(x).T
    y = np.array(y).T
    print(len(x), len(x[0]))
    color = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    fig = plt.figure(figsize=(8,6))	

    plt.xlabel("t [Myr]")
    plt.ylabel("R [pc]")
    for i in range(len(x)):
        plt.plot(x[i], y[i], lw=4, c=color[i])
        plt.scatter(x[i][0], y[i][0], lw=4, c=color[i])
    
    plt.savefig("../figures/IMBHs_inspiral_in_Plummer_potential.pdf")
    plt.show()

def Schwartszchield_radius(mass):
    r = 2 * constants.G * mass/constants.c**2
    return r
    
def main(t_end, M_cl, R_Pl, dt_out, dt_bridge):
    dt_bridge = min(dt_bridge, dt_out)
    #potential = Paczynski_potential(M_cl)
    potential = Plummer_potential(M_cl, R_Pl)
    converter=nbody_system.nbody_to_si(M_cl, R_Pl)
    """
    SMBH = Particle()
    SMBH.mass = M_cl
    SMBH.name = "SMBH"
    SMBH.position = (1,0,0) | units.au
    SMBH.velocity = (0,0,0) | units.kms
    SMBH.radius = Schwartszchield_radius(SMBH.mass)
    SMBH.ncoll = 0
    """

    """
    black_holes = new_plummer_gas_model(5, convert_nbody=converter)
    black_holes.mass = 1000 | units.MSun
    black_holes.radius = Schwartszchield_radius(black_holes.mass)
    black_holes.name = "IMBH"
    black_holes.ncoll = 0
    """

    black_holes = Particles(2)
    black_holes.mass = 50000|units.MSun
    black_holes.name = "IMBH"
    black_holes.position = (10, 0,0) | units.pc
    v = (0.5*constants.G*M_cl/black_holes[0].position.length()).sqrt()
    black_holes.velocity = (0,1,0) *v
    black_holes.radius = Schwartszchield_radius(black_holes.mass)
    black_holes.ncoll = 0

    black_holes[1].position = (5,0,0) | units.pc
    v = (0.1*constants.G*M_cl/black_holes[1].position.length()).sqrt()
    black_holes[1].velocity = (0,1,0) *v
    
    galactic_population = Particles()
    #galactic_population.add_particle(SMBH)
    galactic_population.add_particles(black_holes)

    x, y = integrate_black_holes_in_potential(galactic_population, potential, t_end, dt_out, dt_bridge, converter)

    plot_orbit(x, y)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--dt_out", dest="dt_out", type="float", default = 0.01|units.Myr,
                      help="output_timestep [%default]")
    result.add_option("--dt_bridge", unit=units.yr,
                      dest="dt_bridge", type="float", default = 10000|units.yr,
                      help="bridge timestep [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="M_cl", type="float",default = 1.0e+6 | units.MSun,
                      help="cluster mass [%default]")
    result.add_option("-R", unit= units.parsec,
                      dest="R_Pl", type="float",default = 20.0 | units.parsec,
                      help="Plummer radius [%default]")
    result.add_option("-t", unit= units.Myr,
                      dest="t_end", type="float", default = 10 | units.Myr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

