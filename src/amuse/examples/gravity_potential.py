"""
Gravity potential example from Astrophysical Recipes

Implements a code simulating the galactic center. As the center itself
does not evolve we only need to define the 'get_gravity_at_point' and
'get_potential_at_point'. Note that both functions get arrays of
points.
"""
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.phigrape.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot as plt
from amuse.ic.kingmodel import new_king_model

"""
Implements a code simulating the galactic center. As the center itself
does not evolve we only need to define the 'get_gravity_at_point' and
'get_potential_at_point'. Note that both functions get arrays of
points.
"""

# #BOOKLISTSTART1# #
class GalacticCenterGravityCode(object):
    def __init__(self,R, M, alpha):
        self.radius=R
        self.mass=M
        self.alpha=alpha

    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        m=self.mass*(r/self.radius)**self.alpha  
        fr=constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az

    def circular_velocity(self,r):  
        m=self.mass*(r/self.radius)**self.alpha
        vc=(constants.G*m/r)**0.5
        return vc
# #BOOKLISTSTOP1# #


    def get_potential_at_point(self,eps,x,y,z):
        r=(x**2+y**2+z**2)**0.5
        c=constants.G*self.mass/self.radius**self.alpha    
        phi=c/(self.alpha-1)*(r**(self.alpha-1)-self.radius**(self.alpha-1))
        return phi    
        
# #BOOKLISTSTART3# #
def make_king_model_cluster(nbodycode, N, W0, Mcluster,
                            Rcluster, parameters = []):
      
    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    bodies=new_king_model(N,W0,convert_nbody=converter)

    code=nbodycode(converter)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(bodies)
    return code
# #BOOKLISTSTOP3# #

def plot_cluster(f, x, y, c):
    figure = plt.figure()
    ax = figure.add_subplot(111)
    ax.set_xlabel("X [pc]")
    ax.set_ylabel("Y [pc]")
    #plt.xlim(-60, 60)
    #plt.ylim(-60, 60)
    ax.scatter(x, y, c=c, s=50, lw=0, alpha=0.2)

    save_file = 'Arches.pdf'
    plt.savefig(save_file)
    print('\nSaved figure in file', save_file, '\n')

def evolve_cluster_in_galaxy(N, W0, Rinit, tend, timestep, M, R):

# #BOOKLISTSTART2# #
    Rgal = 1. | units.kpc
    Mgal = 1.6e10 | units.MSun
    alpha = 1.2
    galaxy_code = GalacticCenterGravityCode(Rgal, Mgal, alpha)

    m=galaxy_code.mass*((100|units.pc)/galaxy_code.radius)**galaxy_code.alpha
    cluster_code = make_king_model_cluster(BHTree, N, W0, M, R,
                                           parameters=[("epsilon_squared",
                                                        (0.01 | units.parsec)**2)])
    
    stars = cluster_code.particles.copy()    
    stars.x += Rinit
    stars.vy = 0.8*galaxy_code.circular_velocity(Rinit)
    channel = {"from_framework": stars.new_channel_to(cluster_code.particles),
               "to_framework": cluster_code.particles.new_channel_to(stars)}
    channel["from_framework"].copy_attributes(["x","y","z","vx","vy","vz"])

    plot_cluster(f, stars.x.value_in(units.pc), stars.y.value_in(units.pc), c='r')
    
    system = bridge(verbose=False)
    system.add_system(cluster_code, (galaxy_code,))

    times = quantities.arange(0|units.Myr, tend, timestep)
    xcom = [] | units.pc
    ycom = [] | units.pc
    for i,t in enumerate(times):
        system.evolve_model(t,timestep=timestep)
        channel["to_framework"].copy_attributes(["x","y","z","vx","vy","vz"])
        com = stars.center_of_mass()
        xcom.append(com[0])
        ycom.append(com[1])

    plt.plot(xcom.value_in(units.pc), ycom.value_in(units.pc), lw=1, c='k')
    
    x = system.particles.x.value_in(units.parsec)
    y = system.particles.y.value_in(units.parsec)
    cluster_code.stop()
# #BOOKLISTSTOP2# #

    return x, y


if __name__ == "__main__":
    N=1024
    W0=3
    Rinit=50. | units.parsec
    timestep=0.01 | units.Myr
    #endtime = 2.5 | units.Myr
    endtime = 0.5 | units.Myr
    Mcluster = 5.e4 | units.MSun
    Rcluster = 0.8 | units.parsec
    f = plt.figure(figsize=(8, 8))
    plt.axis('equal')
    plt.scatter([0],[0], marker='+', s=100, lw=1, c='k')
    plt.xlabel('X [pc]')
    plt.ylabel('Y [pc]')
    plt.xlim(-60, 60)
    plt.ylim(-60, 60)

    x, y = evolve_cluster_in_galaxy(N, W0, Rinit, endtime, timestep,
                                    Mcluster, Rcluster)
    plot_cluster(f, x, y, c='b')
    plt.show()
