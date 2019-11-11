/*
 * Interface code
 */

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<sstream>
#include<unistd.h>
#include<unordered_map>
#include<particle_simulator.hpp>
#include<pppt.hpp>
#include"class.hpp"
#include"force.hpp"
#include"init.hpp"
#include"main.hpp"


#ifdef USE_QUAD
using Tree = PS::TreeForForceLong<ForceSoft, EpiSoft, EpjSoft>::QuadrupoleWithSymmetrySearch;
#else
using Tree = PS::TreeForForceLong<ForceSoft, EpiSoft, EpjSoft>::MonopoleWithSymmetrySearch;
#endif
using SystemSoft = PS::ParticleSystem<FpSoft>;

#if defined(USE_PHANTOM)
using CalcForceEpEp = CalcForceEpEpWithLinearCutoffSimd;
using CalcForceEpSp = CalcForceEpSpQuadSimd;
#else
using CalcForceEpEp = CalcForceEpEpWithLinearCutoffNoSimd;
using CalcForceEpSp = CalcForceEpSpQuadNoSimd;
#endif

std::map <int, int> particle_id;

int get_particle_id(int i){
  std::map<int, int>::iterator ii = particle_id.find(i);
  if (ii == particle_id.end()) return -1;
  else return ii->second;
}

void CommBarrier(){
#if defined (USE_BARRIER)
    PS::Comm::barrier();
#endif
}

ParamSet PARAM_SET;

PS::S32 n_glb = 0;
PS::S32 n_loc = 0;
PS::F64 time_begin = 0.0;
PS::F64 time_sys = 0.0;
PS::S32 n_grp_limit = 64;
PS::S32 n_leaf_limit  = 8;
PS::S32 n_smp_ave = 64;
PS::F64 theta = 0.4;

SystemSoft system_soft;

PS::DomainInfo dinfo;

PPPT::ClusterSystem<FpSoft> system_cluster;

// This is causing some trouble
Energy eng_init;
Energy eng_diff_max;

Tree tree_soft;

// ---- AMUSE specific variables

// will be set to 1 the first time we start evolve
bool started_running;
// Counter for all particles. This number must never decrease!
PS::S64 unique_particles = 0;

// ---- end AMUSE specific variables

int initialize_code(){
    started_running = 0;
    PPPT::Initialize();

    //PARAM_SET.time_end_ = 10.0;
    PARAM_SET.time_end_ = 0.125;
    PARAM_SET.eps_    = 0.1;
    PARAM_SET.eps_sq_ = PARAM_SET.eps_ * PARAM_SET.eps_;

    PARAM_SET.theta_    = theta;
    PARAM_SET.dt_soft_  = 1.0 / 256.0;
    
    PARAM_SET.eta_    = 0.1;
    PARAM_SET.eta_s_  = PARAM_SET.eta_ * 0.01;
    PARAM_SET.dt_hard_limit_ = PARAM_SET.dt_soft_ / 8;

    system_soft.initialize();
    system_soft.setAverageTargetNumberOfSampleParticlePerProcess(n_smp_ave);

    dinfo.initialize();    
    return 0;
}

int set_eps2(double epsilon_squared){
    PARAM_SET.eps_sq_ = epsilon_squared;
    PARAM_SET.eps_ = sqrt(epsilon_squared);
    return 0;
}

int get_eps2(double * epsilon_squared){
    *epsilon_squared = PARAM_SET.eps_sq_;
    // eps is not needed
    return 0;
}

int get_time(double * time){
    *time = time_sys;
    return 0;
}

int get_number_of_particles(int * number_of_particles){
    *number_of_particles = n_glb;
    return 0;
}

int get_eta(double * eta){
    *eta = PARAM_SET.eta_; // Hermite timestep criterion
    return 0;
}

int set_eta(double eta){
    PARAM_SET.eta_ = eta;
    PARAM_SET.eta_s_ = PARAM_SET.eta_ * 0.01;
    return 0;
}

int get_theta_for_tree(double * theta_for_tree){
    *theta_for_tree = PARAM_SET.theta_; // Opening angle
    return 0;
}

int set_theta_for_tree(double theta_for_tree){
    PARAM_SET.theta_ = theta_for_tree;
    return 0;
}

int get_radius(int id, double* radius){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    *radius = system_soft[i].r_out;
    return 0;
}

int set_radius(int id, double radius){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    system_soft[i].r_out = radius;
    return 0;
}

int cleanup_code(){
    return 0;
}

int evolve_model(double time){
    started_running = 1;
    Energy eng_init(system_soft);
    std::cerr<<"eng_init:"<<std::endl;
    eng_init.dump(std::cerr);
    
    Energy eng_diff_max;
    PARAM_SET.time_end_ = time;
    PS::S64 n_loop = 0;

    // std::cerr<<"GLOBALS:"<<std::endl;
    // std::cerr<<"n_glb = "<<n_glb<<std::endl;
    // std::cerr<<"n_loc = "<<n_loc<<std::endl;
    // std::cerr<<"time_sys = "<<time_sys<<std::endl;
    // std::cerr<<"n_grp_limit = "<<n_grp_limit<<std::endl;
    // std::cerr<<"n_leaf_limit = "<<n_leaf_limit<<std::endl;
    // std::cerr<<"n_smp_ave = "<<n_smp_ave<<std::endl;
    // std::cerr<<"theta = "<<theta<<std::endl;
    // std::cerr<<std::endl;
    // std::cerr<<"PARAMETERS:"<<std::endl;
    // std::cerr<<"time_end = "<<PARAM_SET.time_end_<<std::endl;
    // std::cerr<<"eps = "<<PARAM_SET.eps_<<std::endl;
    // std::cerr<<"eps_sq = "<<PARAM_SET.eps_sq_<<std::endl;
    // std::cerr<<"theta = "<<PARAM_SET.theta_<<std::endl;
    // std::cerr<<"dt_soft = "<<PARAM_SET.dt_soft_<<std::endl;
    // std::cerr<<"eta = "<<PARAM_SET.eta_<<std::endl;
    // std::cerr<<"eta_s = "<<PARAM_SET.eta_s_<<std::endl;
    // std::cerr<<"dt_hard_limit = "<<PARAM_SET.dt_hard_limit_<<std::endl;
    while(time_sys < PARAM_SET.time_end_){
        std::cerr<<"time_sys:"<<time_sys<<std::endl;
        time_sys += PARAM_SET.dt_soft_;
        std::cerr<<"time_sys ->"<<time_sys<<std::endl;
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"n_loop= "<<n_loop<<std::endl;

        if(n_loop % 4 == 0){
            dinfo.decomposeDomainAll(system_soft);
        }

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 0"<<std::endl;
        
        system_soft.exchangeParticle(dinfo);
        n_loc = system_soft.getNumberOfParticleLocal();
        SetRankToFpOmp(&system_soft[0], n_loc);    

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 1"<<std::endl;
        
        tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSimd(), CalcForceEpSpQuadNoSimd(), system_soft, dinfo);
#if 0
#pragma omp parallel for
        for(auto i=0; i<system_soft.getNumberOfParticleLocal(); i++){
            system_soft[i].pot_tot += system_soft[i].mass / sqrt(system_soft[i].getROut()*system_soft[i].getROut()); // correction for self gravity
        }
#endif

        
        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 2"<<std::endl;
        
        system_cluster.searchClusterLocalAndCorrect(system_soft, tree_soft, CorrectPotAndAcc<FpSoft, EpjSoft>);

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 3"<<std::endl;
        
        KickOmp(system_soft, PARAM_SET.dt_soft_*0.5); // 1st kick


        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 4"<<std::endl;
        
        Energy eng_fin(system_soft);
        Energy eng_diff = (eng_init-eng_fin)/eng_init;
        if(n_loop==0){
            eng_diff_max = eng_diff;
        }
        else{
            if(fabs(eng_diff.getTot()) > fabs(eng_diff_max.getTot())){
                eng_diff_max = eng_diff;
            }
        }
        if(PS::Comm::getRank()==0){
            std::cerr<<"eng_diff:"<<std::endl;
            eng_diff.dump(std::cerr);
            std::cerr<<"eng_diff_max:"<<std::endl;
            eng_diff_max.dump(std::cerr);
        }

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 5"<<std::endl;
        
        KickOmp(system_soft, PARAM_SET.dt_soft_*0.5); // 2nd kick

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 6"<<std::endl;
        
        system_cluster.connectClusterAndIntegrate(system_soft, tree_soft, DriveOneCluster<FpSoft>, DriveMulCluster<FpSoft>, false);

        PS::Comm::barrier();
        if(PS::Comm::getRank()==0) std::cerr<<"check 7"<<std::endl;
        
        n_loop++;
        //if(n_loop == 50) break;
    }

    return 0;
}

int get_position(int id, double* x, double* y, double* z){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    *x = system_soft[i].pos.x;
    *y = system_soft[i].pos.y;
    *z = system_soft[i].pos.z;
    return 0;
}

int get_velocity(int id, double* vx, double* vy, double* vz){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    *vx = system_soft[i].vel.x;
    *vy = system_soft[i].vel.y;
    *vz = system_soft[i].vel.z;
    return 0;
}

int new_particle(int* id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius){
    FpSoft particle;
    PS::S64 pid = unique_particles;
    unique_particles++;
    
    //system_soft.setNumberOfParticleLocal(n_loc);
    //PS::S32 i = unique_particles;

    particle.id = pid;
    particle.mass = mass;
    particle.pos.x = x;
    particle.pos.y = y;
    particle.pos.z = z;
    particle.vel.x = vx;
    particle.vel.y = vy;
    particle.vel.z = vz;
    particle.r_out = radius; // 2./256.;  // PARAM_SET.dt_soft_ * 2.0;  // cutoff radius for ptree-pp
    particle.r_in  = system_soft[i].r_out * 0.1; // inner cutoff radius for ptree-pp
    const PS::F64 r_buf = 1.5 * radius;
    particle.r_search  = system_soft[i].r_out + r_buf;
    //particle.r_out = system_soft[i].r_in = system_soft[i].r_search = 0.0;

    particle_id[pid] = system_soft.getNumberOfParticleGlobal();

    system_soft.addOneParticle(particle);
    *id = pid;
    n_glb++;
    n_loc++;
    //std::cerr<<"n_glb="<<n_glb<<std::endl;
    //std::cerr<<"n_loc="<<n_glb<<std::endl;
    return 0;
}
    
int set_position(int id, double x, double y, double z){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    system_soft[i].pos.x = x;
    system_soft[i].pos.y = y;
    system_soft[i].pos.z = z;
    return 0;
}

int set_velocity(int id, double vx, double vy, double vz){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    system_soft[i].vel.x = vx;
    system_soft[i].vel.y = vy;
    system_soft[i].vel.z = vz;
    return 0;
}

int get_potential(int id, double* potential){
    return -1;    
}

int set_time_step(double time_step){
    PARAM_SET.dt_soft_ = time_step; // timestep for tree code
    return 0;
}

int get_time_step(double* time_step){
    *time_step = PARAM_SET.dt_soft_;
    return 0;
}

int get_begin_time(double* begin_time){
    *begin_time = time_begin;
    return 0;
}

int get_total_mass(double* total_mass){
    return -1;    
}

int set_begin_time(double begin_time){
    if (started_running){
        time_sys = begin_time;
        return 0;
    }
    else {
        return -2; // not allowed to change once we are running
    }
}

int delete_particle(int id){
    PS::S32 i = get_particle_id(id);
    FpSoft p = system_soft[i];
    PS::S32 n_remove = 1;
    PS::S32 *idx[n_remove] = { new PS::S32[n_remove] };
    *idx[0] = p.id;
    system_soft.removeParticle(*idx, n_remove);
    particle_id.erase(id);
    return 0;
}

int commit_particles(){
    dinfo.decomposeDomainAll(system_soft);
    system_soft.exchangeParticle(dinfo);
    n_loc = system_soft.getNumberOfParticleLocal();
    SetRankToFpOmp(&system_soft[0], n_loc);

    PS::Comm::barrier();

    // This bit only needs to be done once
    if(PS::Comm::getRank()==0)
        std::cerr<<"start tree"<<std::endl;
    tree_soft.initialize(n_glb, theta, n_leaf_limit, n_grp_limit);
    tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSimd(), CalcForceEpSpQuadNoSimd(), system_soft, dinfo);

    PS::Comm::barrier();
    if(PS::Comm::getRank()==0)
        std::cerr<<"start cluster search"<<std::endl;
    
    // Also only once
    system_cluster.initialize();
    system_cluster.searchClusterLocalAndCorrect(system_soft, tree_soft, CorrectPotAndAcc<FpSoft, EpjSoft>);

    Energy eng_init(system_soft); // Modified
    if(PS::Comm::getRank()==0){
        std::cerr<<"eng_init:"<<std::endl;
        eng_init.dump(std::cerr);
    }
    
    // Only once
    KickOmp(system_soft, PARAM_SET.dt_soft_*0.5);

    PS::Comm::barrier();
    if(PS::Comm::getRank()==0)
        std::cerr<<"start connect cluster and integrate"<<std::endl;
    
    // Only once
    system_cluster.connectClusterAndIntegrate(system_soft, tree_soft, DriveOneCluster<FpSoft>, DriveMulCluster<FpSoft>, false);

    PS::Comm::barrier();
    if(PS::Comm::getRank()==0)
        std::cerr<<"finish connect cluster and integrate"<<std::endl;

    return 0;
}

int get_acceleration(int id, double* ax, double* ay, double* az){
    return -1;    
}

int get_total_radius(double* total_radius){
    return -1;    
}

int set_acceleration(int id, double ax, double ay, double az){
    return -1;    
}

int commit_parameters(){
    return 0;
}

int synchronize_model(){
    return 0;
}

int get_kinetic_energy(double* kinetic_energy){
    //Energy eng_fin(system_soft);
    //*kinetic_energy = eng_fin.kin; // cannot do this since it is private
    return 0;
}

int recommit_particles(){
    // Do we need to do anything if particles are added, before we use them to
    // evolve further?
    // That goes here.
    return 0;
}

int recommit_parameters(){
    // Do we need to do anything if parameters are changed, before we start the
    // evolve loop?
    // That goes here.
    return 0;
}

int get_potential_energy(double* potential_energy){
    //Energy eng_fin(system_soft);
    //*potential_energy = eng_fin.pot; // cannot do this since it is private
    return 0;
}

int get_index_of_next_particle(int id, int* index_of_next_particle){
    *index_of_next_particle = id+1; // unless it's the last particle
    return 0;
}

int get_center_of_mass_position(double* x, double* y, double* z){
    return -1;    
}

int get_center_of_mass_velocity(double* vx, double* vy, double* vz){
    return -1;    
}

int get_index_of_first_particle(int* id){
    // Unless this particle is deleted?
    *id = 0;
    return 0;
}

int get_mass(int id, double* mass){
    PS::S32 i = get_particle_id(id);
    *mass = system_soft[i].mass;
    return 0;
}

int set_mass(int id, double mass){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    system_soft[i].mass = mass;    
    return 0;
}

int get_state(int id, double* mass, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    *mass = system_soft[i].mass;
    *x = system_soft[i].pos.x;
    *y = system_soft[i].pos.y;
    *z = system_soft[i].pos.z;
    *vx = system_soft[i].vel.x;
    *vy = system_soft[i].vel.y;
    *vz = system_soft[i].vel.z;
    *radius = system_soft[i].r_out;
    return 0;
}

int set_state(int id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius){
    PS::S32 i = get_particle_id(id);
    //PS::S32 i = id;
    system_soft[i].mass = mass;
    system_soft[i].pos.x = x;
    system_soft[i].pos.y = y;
    system_soft[i].pos.z = z;
    system_soft[i].vel.x = vx;
    system_soft[i].vel.y = vy;
    system_soft[i].vel.z = vz;
    system_soft[i].r_out = radius;
    return 0;
}
