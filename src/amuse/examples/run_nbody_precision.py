from amuse.lab import * 
import numpy
import matplotlib.pyplot as plt
import time
import pickle
from amuse.community.hermite_grx.interface import HermiteGRX


def energy_error_of_integrated_Nbody_system(code, particles, end_time, precision):

    gravity = code()
    #if "BHTree" in str(code):
    if code is BHTree:
        print('BHTree')
        gravity.set_time_step(precision | nbody_system.time)
        #gravity.epsilon_squared = 0 | (nbody_system.length)**2
        gravity.parameters.opening_angle = 0.75
        #gravity.parameters.opening_angle = 0.01
        #gravity.parameters.timestep_accuracy_parameter = precision
        #gravity.set_timestep(precision | nbody_system.time)
        #gravity.commit_parameters()
    elif "hermite_grx" in str(code): 
        gravity.parameters.perturbation = '2.5PN_EIH'
        gravity.parameters.integrator = 'Hermite'
        #gravity.parameters.integrator = 'RegularizedHermite'
        #gravity.parameters.integrator = 'SymmetrizedHermite'
        gravity.parameters.light_speed = 100*particles.velocity.max()
        #gravity.parameters.light_speed = 1e+6 | units.length/units.time
        gravity.parameters.dt_param = precision
    else:
        gravity.parameters.timestep_parameter = precision

    gravity.particles.add_particles(particles)
    channel_from_to_framework = gravity.particles.new_channel_to(particles)

    E0 = gravity.particles.potential_energy(G=nbody_system.G) 
    E0 += gravity.particles.kinetic_energy()
    t_cpu = time.time()
    gravity.evolve_model(end_time)
    t_cpu = time.time()-t_cpu 
    channel_from_to_framework.copy()
    Et = gravity.particles.potential_energy(G=nbody_system.G) \
                            + gravity.particles.kinetic_energy()
    gravity.stop()

    de = -abs(Et-E0)/E0
    return de, t_cpu|units.s

def get_dE(code, precision, t_end, filename):
    numpy.random.seed(31415)
    particles = new_plummer_model(1000)
    
    dE = []
    dt = [] | units.s
    for pri in precision:
        dEi, t_cpu = energy_error_of_integrated_Nbody_system(code,
                                                             particles,
                                                             t_end, pri)
        dE.append(abs(dEi))
        dt.append(t_cpu)
        print("integrated", code, "with precision=",
              pri, "dE/E=", dEi, "t_cpu=", t_cpu.in_(units.s))
        pickle.dump([precision, dE, dt] ,
                    open(filename, "wb"))
        
    return dE, dt
    
if __name__ in ('__main__','__plot__'):

    precision = 10.**numpy.linspace(0.9, -3., 13)
    #precision = 10.**numpy.linspace(1., -2., 3)

    t_end = 1.0| nbody_system.time
    cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    figure = plt.figure(figsize=(6,6))
    ax = figure.add_subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')

    """
    print('ph4')
    code = ph4
    dE, dt = get_dE(code, precision, t_end, "ph4_precision.pkl")
    plt.scatter(precision, dE, c=cols[0], lw=0, s=100, marker='o')
    #pickle.dump([precision, dE, dt] , open( "ph4_precision.pkl", "wb" ) )
    """

    code = BHTree
    dE, dt = get_dE(code, precision, t_end, "BHTree_precision.pkl")
    plt.scatter(precision, dE, c=cols[1], lw=0, s=100, marker='^')
    
    """
    print('Hermite-GRX')
    code = HermiteGRX
    dE, dt = get_dE(code, precision, 0.01*t_end, "Hermitegrx_precision.pkl")
    plt.scatter(precision, dE, c=cols[0], lw=0, s=100, marker='o')
    import pickle
    #pickle.dump([dE, dt] , open( "Hermitegrx_precision.pkl", "wb" ) )
    """

    
