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
//Energy eng_init(system_soft);

//Energy eng_diff_max;
PS::S64 n_loop = 0;

Tree tree_soft;

int initialize_code(){
    PPPT::Initialize();

    //PARAM_SET.time_end_ = 10.0;
    PARAM_SET.time_end_ = 0.125;
    PARAM_SET.eps_    = 0.01;
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
    *eta = PARAM_SET.eta_;
    return 0;
}

int set_eta(double eta){
    PARAM_SET.eta_ = eta;
    return 0;
}

int get_theta_for_tree(double * theta_for_tree){
    *theta_for_tree = PARAM_SET.theta_;
    return 0;
}

int set_theta_for_tree(double theta_for_tree){
    PARAM_SET.theta_ = theta_for_tree;
    return 0;
}

int get_radius(int, double* radius){
    return 0;
}

int set_radius(int, double){
    return 0;
}

int cleanup_code(){
    return 0;
}

int evolve_model(double time){
    Energy eng_init(system_soft);
    
    Energy eng_diff_max;
    PARAM_SET.time_end_ = time;
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
    PS::S32 i = id;
    *x = system_soft[i].pos.x;
    *y = system_soft[i].pos.y;
    *z = system_soft[i].pos.z;
    return 0;
}

int get_velocity(int id, double* vx, double* vy, double* vz){
    PS::S32 i = id;
    *vx = system_soft[i].vel.x;
    *vy = system_soft[i].vel.y;
    *vz = system_soft[i].vel.z;
    return 0;
}

int new_particle(int* id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius){
    n_glb += 1;
    PS::S32 my_rank = PS::Comm::getRank();
    PS::S32 n_proc = PS::Comm::getNumberOfProc();
    if(my_rank == 0){
        n_loc = n_glb;
    }
    else{
        n_loc = 0;
    }
    system_soft.setNumberOfParticleLocal(n_loc);
    PS::S32 i = n_glb;

    system_soft[i].mass = mass;
    system_soft[i].pos.x = x;
    system_soft[i].pos.y = y;
    system_soft[i].pos.z = z;
    system_soft[i].vel.x = vx;
    system_soft[i].vel.y = vy;
    system_soft[i].vel.z = vz;
    system_soft[i].id = i;
    *id = i;
    return 0;
}
    
int set_position(int id, double x, double y, double z){
    PS::S32 i = id;
    system_soft[i].pos.x = x;
    system_soft[i].pos.y = y;
    system_soft[i].pos.z = z;
    return 0;
}

int set_velocity(int id, double vx, double vy, double vz){
    PS::S32 i = id;
    system_soft[i].vel.x = vx;
    system_soft[i].vel.y = vy;
    system_soft[i].vel.z = vz;
    return 0;
}

int get_potential(int id, double* potential){
    return -1;    
}

int set_time_step(double time_step){
    PARAM_SET.dt_soft_ = time_step;
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
    return -1;    
}

int delete_particle(int id){
    return -1;    
}

int commit_particles(){
    dinfo.decomposeDomainAll(system_soft);
    system_soft.exchangeParticle(dinfo);
    n_loc = system_soft.getNumberOfParticleLocal();
    SetRankToFpOmp(&system_soft[0], n_loc);
    for(auto i=0; i<n_loc; i++){
        system_soft[i].r_out = PARAM_SET.dt_soft_ * 2.0;
        system_soft[i].r_in  = system_soft[i].r_out * 0.1;
        const PS::F64 r_buf = 3.0*PARAM_SET.dt_soft_;
        system_soft[i].r_search  = system_soft[i].r_out + r_buf;
        //system_soft[i].r_out = system_soft[i].r_in = system_soft[i].r_search = 0.0;
    }

    PS::Comm::barrier();

        if(PS::Comm::getRank()==0)
        std::cerr<<"start tree"<<std::endl;
    tree_soft.initialize(n_glb, theta, n_leaf_limit, n_grp_limit);
    tree_soft.calcForceAllAndWriteBack(CalcForceEpEpWithLinearCutoffNoSimd(), CalcForceEpSpQuadNoSimd(), system_soft, dinfo);

#if 0
#pragma omp parallel for
    for(auto i=0; i<system_soft.getNumberOfParticleLocal(); i++){
        system_soft[i].pot_tot += system_soft[i].mass / sqrt(system_soft[i].getROut()*system_soft[i].getROut()); // correction for self gravity
    }
#endif
    
#if 0
    PS::ReallocatableArray<ForceSoft> force_dir(n_loc, n_loc, PS::MemoryAllocMode::Pool);
    for(auto i=0; i<n_loc; i++)
        force_dir[i].clear();
    tree_soft.calcForceDirectParallel
        (CalcForceNewtonNoSimd(), &force_dir[0], system_soft, dinfo);
    Energy eng_dir;
    eng_dir.set(system_soft, force_dir);
    if(PS::Comm::getRank()==0){
        std::cerr<<"eng_dir:"<<std::endl;
        eng_dir.dump(std::cerr);
    }

    Energy eng_tmp(system_soft);
    if(PS::Comm::getRank()==0){
        std::cerr<<"eng_tmp:"<<std::endl;
        eng_tmp.dump(std::cerr);
    }
#endif


    PS::Comm::barrier();
    if(PS::Comm::getRank()==0)
        std::cerr<<"start cluster search"<<std::endl;
    
    system_cluster.initialize();
    system_cluster.searchClusterLocalAndCorrect(system_soft, tree_soft, CorrectPotAndAcc<FpSoft, EpjSoft>);

    //FIXME if(PS::Comm::getRank()==0){
    //FIXME     std::cerr<<"eng_init:"<<std::endl;
    //FIXME     eng_init.dump(std::cerr);
    //FIXME }
    
    KickOmp(system_soft, PARAM_SET.dt_soft_*0.5);

    PS::Comm::barrier();
    if(PS::Comm::getRank()==0)
        std::cerr<<"start connect cluster and integrate"<<std::endl;
    
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
    return -1;    
}

int recommit_particles(){
    return 0;
}

int recommit_parameters(){
    return 0;
}

int get_potential_energy(double* potential_energy){
    return -1;    
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
    *id = 0;
    return 0;
}

int get_mass(int id, double* mass){
    PS::S32 i = id;
    *mass = system_soft[i].mass;
    return 0;
}

int set_mass(int id, double mass){
    PS::S32 i = id;
    system_soft[i].mass = mass;    
    return 0;
}

int get_state(int id, double* mass, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius){
    PS::S32 i = id;
    *mass = system_soft[i].mass;
    *x = system_soft[i].pos.x;
    *y = system_soft[i].pos.y;
    *z = system_soft[i].pos.z;
    *vx = system_soft[i].vel.x;
    *vy = system_soft[i].vel.y;
    *vz = system_soft[i].vel.z;
    return 0;
}

int set_state(int id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius){
    PS::S32 i = id;
    system_soft[i].mass = mass;
    system_soft[i].pos.x = x;
    system_soft[i].pos.y = y;
    system_soft[i].pos.z = z;
    system_soft[i].vel.x = vx;
    system_soft[i].vel.y = vy;
    system_soft[i].vel.z = vz;
    return 0;
}
