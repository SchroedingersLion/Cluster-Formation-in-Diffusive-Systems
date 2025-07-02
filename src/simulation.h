#ifndef SIMULATION_H
#define SIMULATION_H

#include "model.h"
#include "measurement.h"
#include "coordinate.h"

#include <omp.h>
#include <random>
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>

#define _USE_MATH_DEFINES

/*
Holds the simulation class definition. This class applies an integrator to propagate the positions, velocities and forces
of the system stored in the passed model. The class can easily be extended by additional integrators. While the details of the
pairwise forces between particles are specified in the model class, this class implements a multithreaded overall force routine
that iterates through all pairs in the system and calls the pairwise force function of the model instance.
While simulating a system, the simulation class regularly accesses the measurement class to store current observable values.

All member functions are implemented inline for simplicity.
*/


// ###################### SIMULATION CLASS DEFINITION ##################################################
template <size_t DIMENSION>
class simulation {

    public: 

        measurement<DIMENSION> meas;  // Needs to be public because main() calls print function of the measurement class.

        // CONSTRUCTOR.
        simulation<DIMENSION> (IPS_model<DIMENSION>& model, 
                    const double stepsize, 
                    const double beta, 
                    const double gamma, 
                    const int N_iter,
                    const int N_meas, 
                    const int threads, 
                    const std:: string integrator, 
                    const std:: string init_conf, 
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
                    init_conf {init_conf}, 
                    seed {seed},
                    meas {measurement (model, N_meas, N_iter, stepsize, trajectory)}
            {

                std:: cout << "\nCreating simulation with integrator " << integrator << ",\n"
                           << "stepsize " << stepsize << ",\n"
                           << "gamma " << gamma << ",\n"
                           << "beta " << beta << ",\n"
                           << N_iter << " iterations.\n"
                           << "Using " << THREADS << " threads for force calculations. \n"
                           << "Randomseed " << seed << ".\n"
                           << "Taking a measurement any " << N_meas << " steps.\n" << std:: endl;

                // Specify integrator.
                if (integrator=="baoab") integrator_step = &simulation<DIMENSION>:: BAOAB_step;
                else if (integrator=="ubu") integrator_step = &simulation:: UBU_step;
                else throw std:: invalid_argument( "Invalid integrator argument. Allowed are 'baoab' and 'ubu'." );

                // Help vector needed for parallel force computation.
                forces_for_all_tasks.resize(THREADS);
                std:: fill(forces_for_all_tasks.begin(), forces_for_all_tasks.end(), std:: vector <coordinate<DIMENSION>> (model.N_particles, coordinate<DIMENSION>()));

                twister.seed(seed);       // Seed RNG for simulation.

            }; 

        void run();   // Run simulation.



    private:
        
        std:: mt19937 twister;                    // RNG used in the integrators.
        std:: normal_distribution<> normal{0,1};  // Distribution used in the integrators.
        IPS_model<DIMENSION>& model;              // Stores the particle system.
        
        // Integrator parameters. 
        const double stepsize;          
        const double beta;
        const double gamma;
        const int N_iter;
        const int N_meas;
        const int seed;
        const std:: string integrator;  // Which integrator to use.

        /* 
        How to initialize positions. Either "uniform" or name of .csv file storing 
        two columns (delimiter " ") of particle positions.
        */
        const std:: string init_conf;   

        const int THREADS;  // Number of threads used to compute forces.
        std:: vector <std:: vector <coordinate<DIMENSION>>> forces_for_all_tasks;   // One force vector per thread.
        


        // Integrator functions.
        void (simulation<DIMENSION>::* integrator_step)();   // Point to integrator to be used.
        void BAOAB_step();
        void UBU_step();
        void A_step(const double h);
        void B_step(const double h);
        void O_step(const double h);
        void U_step(const double h);
        void compute_force_par();           // Multithreaded force routine.
        void apply_periodic_boundaries();

        // Functions for system initialization.
        std::vector<std::vector<double>>  read_from_file(const std:: string& filename);
        void set_initial_position(const std:: string& init_conf);
        void set_initial_position(const int seed);
        void set_initial_velocities();

};
// ##################### END OF CLASS DEFINITION ##############################################


// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################
template <size_t DIMENSION>
inline void simulation<DIMENSION>:: run(){
    /*
    Main function to run the simulation. Initializes the system, seeds the RNG, and executes the
    integrator iterations (taking observable measurements in between).
    */

    // Set positions.
    if (!init_conf.empty() && init_conf != "uniform") set_initial_position(init_conf); 
    else set_initial_position(seed);

    set_initial_velocities(); // Set velocities.

    compute_force_par();      // Set forces.


    // Main loop.
    std:: normal_distribution<> normal{0,1};
    const double h_half = 0.5 * stepsize;
    
    std:: cout << "\nRunning simulation...\n";

    auto t1 = std:: chrono:: high_resolution_clock:: now();  // Measure runtime.
    
    for (int i=0; i<=N_iter; ++i){
        
        if (i % N_meas == 0) meas.take_measurement();

        (this->*integrator_step)();

        if (i%1000==0) std:: cout << "Iteration "<<i<< " done!" << std:: endl;

    }

    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: UBU_step(){

    const double h_half = stepsize/2;

    U_step(h_half);
    apply_periodic_boundaries();
    compute_force_par();
    B_step(stepsize);
    U_step(h_half);

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: BAOAB_step(){

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


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: A_step(const double h){

    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<DIMENSION; ++dim){
            model.positions[i][dim] += h*model.velocities[i][dim];
        }

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: B_step(const double h){
    
    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<DIMENSION; ++dim){
            model.velocities[i][dim] += h*model.forces[i][dim];
        }

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: O_step(const double h){

    
    const double a = exp(-gamma*h);
    const double pref = sqrt(1/beta *(1-a*a));
    /* 
    Note: Recreating these constants is somewhat inefficient, but they can't be precomputed 
    if O_step should be callable for different step-widths h, which is the case here.
    Passing h as an argument allows for the usage of the same O_step code for integrators 
    other than BAOAB (e.g., OBABO) to be added in the future. 
    Need a cleaner way to handle these constants without using flexibility.
    */

    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<DIMENSION; ++dim){
            model.velocities[i][dim] = a*model.velocities[i][dim] + pref*normal(twister); 
        }

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: U_step(const double h){

    // Prefactors used in the U step.
    const double pref_U1 {exp(-gamma*h)};
    const double pref_U2 {(1-pref_U1)/gamma};
    const double pref_U3 {sqrt(2/(beta*gamma))};
    const double pref_U4 {sqrt(2*gamma/beta)};
    const double pref_Z1 {sqrt(h)};
    const double pref_Z2 {sqrt( (1-pref_U1*pref_U1)/(2*gamma) )};
    const double pref_Z3 {sqrt( 2*(1-pref_U1)/(gamma*h*(1+pref_U1)) )};
    const double pref_Z4 {sqrt( 1-2*(1-pref_U1)/(gamma*h*(1+pref_U1)) )};
    /* 
    Note: Recreating these constants is inefficient (see note in O_step function above).
    */

    const double pref_Z2_total1 {pref_Z2 * pref_Z3};
    const double pref_Z2_total2 {pref_Z2 * pref_Z4};


    // Compute new positions/velocities.

    double xi1, xi2, Z1, Z2; 

    for (int i=0; i<model.N_particles; ++i){

        for (size_t dim = 0; dim<DIMENSION; ++dim){

            xi1 = normal(twister);
            xi2 = normal(twister);

            Z1 = pref_Z1 * xi1;

            Z2 = pref_Z2_total1 * xi1 + pref_Z2_total2 * xi2;

            model.positions[i][dim] += pref_U2 * model.velocities[i][dim] + pref_U3 * (Z1 - Z2);

            model.velocities[i][dim] = pref_U1 * model.velocities[i][dim] + pref_U4 * Z2;
        
        }

    }

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: compute_force_par(){
    /*
    Fills the forces vector of the model instance by calling the pairwise interaction routine
    specified in model. To do this, this function needs to iterate over all particle pairs, which is
    an O(N^2) operation. The iterations are distributed across different threads managed via OpenMP.
    */

    // Initialize force vector on each task (will have shape (N,DIMENSION)).
    for (auto& forces_for_specific_task : forces_for_all_tasks)
        std:: fill(forces_for_specific_task.begin(), forces_for_specific_task.end(), coordinate<DIMENSION>());

    // OMP environment.
    #pragma omp parallel num_threads(THREADS)
    {
    
        coordinate <DIMENSION> distance;    // Store distance (per dimension) between two particles.
        coordinate <DIMENSION> force_ij;    // Store force (per dimension) between those two particles.

        // Dynamically distribute iterations over particle paris across OMP threads.
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < model.N_particles; ++i) {
            for (int j = i + 1; j < model.N_particles; ++j) {
                
                // Fill distance and force_ij.
                model.get_distances_ij(i, j, distance);             
                (model.*(model.get_force_ij))(distance, force_ij); 

                // Access thread specific force vector.
                std:: vector <coordinate<DIMENSION>>& forces_for_this_task = forces_for_all_tasks[omp_get_thread_num()];
                
                // Add forces to thread specific force vector.
                for (size_t dim = 0; dim<DIMENSION; ++dim){
                    forces_for_this_task[i][dim] += force_ij[dim];
                    forces_for_this_task[j][dim] += -force_ij[dim];
                }

            }
        }
    }


    std:: fill(model.forces.begin(), model.forces.end(), coordinate<DIMENSION>());   // Set forces in model to 0.

    // Sum all of the task-specific forces into the output parameter.
    for (auto const& forces_for_specific_task : forces_for_all_tasks)
        for (int i = 0; i < model.N_particles; ++i)
            for (size_t dim = 0; dim<DIMENSION; ++dim){
                model.forces[i][dim] += forces_for_specific_task[i][dim];
            }

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: apply_periodic_boundaries(){
    
    double pos;
    const double L {model.L};
    const double L_half {0.5*L};

    for (int i=0; i<model.N_particles; ++i)
        for (size_t dim = 0; dim<DIMENSION; ++dim){
            pos = model.positions[i][dim];
            model.positions[i][dim] = pos>L_half ? pos-L : (pos<-L_half ? pos + L : pos);
        }

}


template <size_t DIMENSION>
inline std::vector<std::vector<double>> simulation<DIMENSION>:: read_from_file(const std:: string& filename){
    /* 
    Read .csv file filename containing unknown number of rows and columns. 
    Each row MUST consist of numbers separated by a single space symbol " ". 
    */

    std::vector<std::vector<double>> data;  // Store data read from file.
    
    std::ifstream file(filename);   // Open file.
    
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        
        std::vector<double> row;    // Store single row from file.
        std::string cell;           // Store single element from row.

        std::stringstream ss(line); // Allow decomposition of one row into individual elements obtainable via getline().

        while (std::getline(ss, cell, ' ')) {
            try {
                row.push_back(std::stod(cell));
            } 
            catch (...) {
                throw std::runtime_error(std::string("Invalid number in CSV: ") + cell);
            }
        }

        // Prevent adding of empty lines.
        if (!row.empty()) {
            data.push_back(std::move(row));   
        }
    }

    return data;
}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: set_initial_position(const std:: string& init_conf){
    // Set initial positions as read from file init_conf.
    
    std:: cout << "Opening file " << init_conf << " to read initial configuration..." << std:: endl;
    std::vector<std::vector<double>> data =  read_from_file(init_conf);
    
    // Check whether data has consistend rows and columns.
    const size_t N_rows {data.size()};
    if (N_rows != model.N_particles) throw std::runtime_error( std:: string("Passed number of particles via --N_particles "
                                                               "does not coincide with number of rows in file ") + init_conf + ".");
    for (auto row : data){
        if (row.size() != DIMENSION) throw std::runtime_error(std:: string("Number of columns in (at least one row) of file ") + init_conf +
                                                              " does not coincide with system dimension passed via --dimension.");
    }

    double pos;
    for (int i=0; i<model.N_particles; ++i)
        for (int dim=0; dim<DIMENSION; ++dim){

            pos = data[i][dim];
            // Check whether read coordinate lies in the simulation box.
            if (pos > model.L/2. || pos < -model.L/2.) throw std::runtime_error(std:: string("One of the coordinates from file ") + init_conf +
                                                                         " does exceed the box boundaries, whose side has length [-L,L] "
                                                                         " where 2L was given via flag --boxlength (default 10, so L=5).");

            model.init_positions[i][dim] = data[i][dim];
        
        }

    model.positions = model.init_positions;

    std:: cout << "Successfully read configuration." << std:: endl;

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: set_initial_position(const int seed){
    // Uniform initialization.

    std:: uniform_real_distribution<double> box_uniform(-model.L/2., model.L/2.);

    for (int i=0; i<model.N_particles; ++i){
        for (size_t dim=0; dim<DIMENSION; ++dim) model.init_positions[i][dim] = box_uniform(twister);
    }

    model.positions = model.init_positions;

}


template <size_t DIMENSION>
inline void simulation<DIMENSION>:: set_initial_velocities(){
    // Gaussian initialization.

    std:: normal_distribution<> normal{0, sqrt(1/beta)};
    
    for (int i=0; i<model.N_particles; ++i){
        for (size_t dim=0; dim<DIMENSION; ++dim) model.velocities[i][dim] = normal(twister);
    }

}

// ########### END OF MEMBER DEFINITIONS ##############################################

#endif // SIMULATION_H



