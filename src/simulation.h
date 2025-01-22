#ifndef SIMULATION_H
#define SIMULATION_H

#define _USE_MATH_DEFINES

#include "model.h"
#include "measurement.h"

#include <omp.h>
#include <random>
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>


// ###################### SIMULATION CLASS DEFINITION ##################################################

class simulation {

    public: 

        measurement meas;  // Needs to be public because main calls print function.

        // CONSTRUCTOR.
        simulation (IPS_model& model, 
                    const double stepsize, 
                    const double beta, 
                    const double gamma, 
                    const int N_iter,
                    const int N_meas, 
                    const int threads, 
                    const std:: string integrator, 
                    const std:: string init_mode, 
                    const int seed,
                    const bool trajectory)
                  : model {model}, 
                    stepsize {stepsize}, 
                    beta {beta}, 
                    gamma {gamma}, 
                    N_iter {N_iter},
                    N_meas {N_meas}, 
                    THREADS {threads}, 
                    integrator {integrator}, 
                    init_mode {init_mode}, 
                    seed {seed},
                    meas {measurement (model, N_meas, N_iter, stepsize, trajectory)}
            {

                std:: cout << "\nCreating simulation with integrator " << integrator << ",\n"
                           << "stepsize " << stepsize << ",\n"
                           << "gamma " << gamma << ",\n"
                           << "beta " << beta << ",\n"
                           << N_iter << " iterations.\n"
                           << "Using " << THREADS << " threads for force calculations. \n"
                           << "Mode of initialization: " << init_mode << ".\n"
                           << "Randomseed " << seed << ".\n"
                           << "Taking a measurement any " << N_meas << " steps.\n" << std:: endl;

                // Specify integrator.
                if (integrator=="BAOAB") integrator_step = &simulation:: BAOAB_step;
                // else if (integrator=="UBU") integrator_step = &simulation:: UBU_step;
                else throw std:: invalid_argument( "Invalid integrator argument. Allowed are 'BAOAB' and 'UBU'." );

                // Check initialization mode.
                if (init_mode!="uniform" && init_mode!="grid") throw std:: invalid_argument("Invalid mode of initialization. Allowed are 'uniform' or 'grid'.");

                // Help vector needed for parallel force computation.
                forces_for_all_tasks.resize(THREADS);
                std:: fill(forces_for_all_tasks.begin(), forces_for_all_tasks.end(), std:: vector <double> (model.N_particles*model.dimension));

            }; 

        void run();



    private:
        
        std:: mt19937 twister;
        IPS_model& model;
        const double stepsize;
        const double beta;
        const double gamma;
        const int N_iter;
        const int N_meas;
        const int THREADS;
        std:: vector <std:: vector <double>> forces_for_all_tasks;
        const std:: string integrator;  
        const std:: string init_mode;
        const int seed;
        void (simulation::* integrator_step)();
        void compute_force_par();
        void A_step(const double h);
        void B_step(const double h);
        void O_step(const double h);
        // void U_step(const double h);
        void apply_periodic_boundaries();
        void BAOAB_step();
        // void UBU_step();
        // void set_initial_position();
        void set_initial_position(const int seed);
        void set_initial_velocities();


};
// ##################### END OF CLASS DEFINITION ##############################################


// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################
// inline void simulation:: set_initial_position(){
//     // Square lattice initialization.

//     // Obtain number of particles per box dimension (for even spacing).
//     int Nx = model.N_particles;
    
//     const double L_mod = model.L-0.00001; // To ensure no particle is placed on the edge.

//     const double dx = 2*L_mod/Nx; // Spacing between particles.

//     // Place particles.
//     for(int ix=0; ix<Nx; ++ix){

//         model.init_positions[ix].x = -L_mod + ix*dx;

//     }

//     model.positions = model.init_positions;

// }



inline void simulation:: set_initial_position(const int seed){
    // Uniform initialization.

    twister.seed(seed);
    std:: uniform_real_distribution<double> box_uniform(-model.L, model.L);
    for (int i=0; i<model.N_particles; ++i){
        for (size_t dim=0; dim<model.dimension; ++dim) model.init_positions[i*model.dimension + dim] = box_uniform(twister);
    }

    model.positions = model.init_positions;

}



inline void simulation:: set_initial_velocities(){
    // Gaussian initialization.

    twister.seed(seed);
    std:: normal_distribution<> normal{0, sqrt(1/beta)};
    
    for (int i=0; i<model.N_particles; ++i){
        for (size_t dim=0; dim<model.dimension; ++dim) model.velocities[i*model.dimension + dim] = normal(twister);
    }

}



inline void simulation:: compute_force_par()
{

    #pragma omp parallel num_threads(THREADS)
    {
        thread_local std::vector<double> distance(model.dimension);
        thread_local std::vector<double> force_ij(model.dimension);

        #pragma omp for schedule(dynamic)
        for (int i = 0; i < model.N_particles; ++i) {
            for (int j = i + 1; j < model.N_particles; ++j) {
                
                model.get_distances_ij(i, j, distance);
                (model.*(model.get_force_ij))(distance, force_ij);

                std:: vector <double>& forces_for_this_task = forces_for_all_tasks[omp_get_thread_num()];
                
                for (size_t dim = 0; dim<model.dimension; ++dim){
                    forces_for_this_task[i*model.dimension + dim] += force_ij[dim];
                    forces_for_this_task[j*model.dimension + dim] += -force_ij[dim];
                }

            }
        }
    }

        // Sum all of the task-specific forces into the output parameter.
    std:: fill(model.forces.begin(), model.forces.end(), 0.0);
    for (auto const& forces_for_specific_task : forces_for_all_tasks)
        for (int i = 0; i < model.N_particles; ++i)
            for (size_t dim = 0; dim<model.dimension; ++dim){
                model.forces[i*model.dimension + dim] += forces_for_specific_task[i*model.dimension + dim];
            }

}



inline void simulation:: A_step(const double h){

    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<model.dimension; ++dim){
            model.positions[i*model.dimension + dim] += h*model.velocities[i*model.dimension + dim];
        }

}



inline void simulation:: B_step(const double h){
    
    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<model.dimension; ++dim){
            model.velocities[i*model.dimension + dim] += h*model.forces[i*model.dimension + dim];
        }

}



inline void simulation:: O_step(const double h){

    std:: normal_distribution<> normal{0,1};
    const double a = exp(-gamma*h);
    const double pref = sqrt(1/beta *(1-a*a));

    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<model.dimension; ++dim){
            model.velocities[i*model.dimension + dim] = a*model.velocities[i*model.dimension + dim] + pref*normal(twister); 
        }

}



// inline void simulation:: U_step(const double h){

// 	std:: normal_distribution<> normal{0,1};

//     // Prefactors used in the U step.
//     const double pref_U1 {exp(-gamma*h)};
//     const double pref_U2 {(1-pref_U1)/gamma};
//     const double pref_U3 {sqrt(2/(beta*gamma))};
//     const double pref_U4 {sqrt(2*gamma/beta)};
//     const double pref_Z1 {sqrt(h)};
//     const double pref_Z2 {sqrt( (1-pref_U1*pref_U1)/(2*gamma) )};
//     const double pref_Z3 {sqrt( 2*(1-pref_U1)/(gamma*h*(1+pref_U1)) )};
//     const double pref_Z4 {sqrt( 1-2*(1-pref_U1)/(gamma*h*(1+pref_U1)) )};

//     const double pref_Z2_total1 {pref_Z2 * pref_Z3};
//     const double pref_Z2_total2 {pref_Z2 * pref_Z4};


//     // Compute new positions/velocities.
//     coordinate Z1, Z2;
//     coordinate xi1, xi2; 

//     for (int i=0; i<model.N_particles; ++i){

//         xi1.x = normal(twister);
//         xi2.x = normal(twister);

//         Z1.x = pref_Z1 * xi1.x;

//         Z2.x = pref_Z2_total1 * xi1.x + pref_Z2_total2 * xi2.x;

//         model.positions[i].x += pref_U2 * model.velocities[i].x + pref_U3 * (Z1.x - Z2.x);

//         model.velocities[i].x = pref_U1 * model.velocities[i].x + pref_U4 * Z2.x;

//     }

// }



inline void simulation:: apply_periodic_boundaries(){
    
    double pos;
    const double L {model.L};
    const double two_L {2*L};

    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<model.dimension; ++dim){
            pos = model.positions[i*model.dimension + dim];
            model.positions[i*model.dimension + dim] = pos>L ? pos-two_L : (pos<-L ? pos + two_L : pos);
        }
}



// inline void simulation:: UBU_step(){

//     const double h_half = stepsize/2;

//     U_step(h_half);
//     apply_periodic_boundaries();
//     compute_force_par();
//     B_step(stepsize);
//     U_step(h_half);

// }



inline void simulation:: BAOAB_step(){

    const double h_half = stepsize/2;

    B_step(h_half);
    A_step(h_half);
    apply_periodic_boundaries();
    O_step(stepsize);
    A_step(h_half);
    apply_periodic_boundaries();
    compute_force_par();
    B_step(h_half);

}



inline void simulation:: run(){

    std:: cout << "\nRunning simulation...\n";
    //  std::cout<< "Simulation at T/Tcrit="<<T_Tcrit_gauss<<" with gamma="<<gamma<<std::endl;

    // Prepare simulation.

    // Set positions.
    if (init_mode == "uniform") set_initial_position(seed);
    // else if (init_mode == "grid") set_initial_position();

    const std:: vector <std:: vector <double> > init_positions {model.positions};  // Needed for MSD computation.

    // std:: fill(model.velocities.begin(), model.velocities.end(), coordinate{0,0});  // Reset velocities.
    set_initial_velocities();

    // Seed RNG for simulation.
    twister.seed(seed);

    // Set forces.
    compute_force_par();


    // Main loop.
    std:: normal_distribution<> normal{0,1};
    const double h_half = 0.5 * stepsize;

    auto t1 = std:: chrono:: high_resolution_clock:: now();
    for (int i=0; i<=N_iter; ++i){
        
        if (i % N_meas == 0) meas.take_measurement();

        (this->*integrator_step)();

        if (i%1000==0) std:: cout << "Iteration "<<i<< " done!" << std:: endl;

    }

    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";

}

// ########### END OF MEMBER DEFINITIONS ##############################################

#endif // SIMULATION_H



