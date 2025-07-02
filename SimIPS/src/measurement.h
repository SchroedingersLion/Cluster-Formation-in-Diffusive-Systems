#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "model.h"
#include "coordinate.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>

#define _USE_MATH_DEFINES

/*
Holds the measurement class which specifies the observables that are collected during simulation (including their computation),
and provides the print function used to print observable time series and (if required) the whole trajectory to files.

New observables can easily be added. For this, make the new observable known to the constructor analogous to the existing ones,
and implement its expression from the system's configuration in a new member function. 
Code regions which the user needs to modify for this aim are marked by 
\\ !!!!!!!

\\ !!!!!!! 
environments.

All member functions are inline for simplicity.
*/


// ###################### MEASUREMENT CLASS DEFINITION ##################################################
template <size_t DIMENSION>
class measurement {

    public:

        // CONSTRUCTOR.
        measurement<DIMENSION>( const IPS_model<DIMENSION>& model, 
                                int N_meas, 
                                const int N_iter, 
                                const double stepsize, 
                                const bool trajectory)
            : model {model}, N_meas {N_meas}, N_iter {N_iter}, stepsize {stepsize}, trajectory {trajectory}
            {
                
                // !!!!!!!!!! ENTER THE NUMBER OF OBSERVABLES TO COLLECT !!!!!!!!!!!!!!!!!!!!!
                no_observables = 3;   // Modify if needed.
                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                observables.resize(no_observables);
                results.resize(no_observables);
                col_names.resize(no_observables);

                // !!!!!!!!!! ENTER NAMES OF OBSERVABLES (WILL BE HEADER OF OUTPUTFILE) !!!!!!!
                col_names[0] = "COM";
                col_names[1] = "MSD";
                col_names[2] = "Tkin";

                // Add new entry of needed.

                // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                

                // Number of measurements that will be taken during simulation based on N_iter and N_meas.
                int no_of_measurements {N_iter / N_meas + 1};
                

                // Resize vectors.
                for (auto& observable_vector : results) observable_vector.resize(no_of_measurements);
                times.resize(no_of_measurements);
                if (trajectory) trajectory_buffer.resize(no_of_measurements);

            
            };

        // Default constructor.
        measurement<DIMENSION>(){};


        void take_measurement(){

            /* !!!!!!!!!!!!! COMPUTE CURRENT OBSERVABLE VALUES FROM PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               The number of entries in vector "observables" must correspond to member variable 
               "no_observables" set by the user in the constructor above. */          
            
            observables[0] = get_center_of_mass_distance();
            observables[1] = get_msd();
            observables[2] = get_Tkin();

            // Add entry if needed.

            // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            add_to_results();  // Add taken observables to results array and add new time value.

            return;

        };


        void print_results(const std:: string outputname);    // Print results array to .csv.


    private:
        const IPS_model<DIMENSION>& model;             // Model to take measurements on.
        int no_observables;                            // Number of observables to be taken.
        std:: vector <float> observables;              // Vector of size (no_observables) storing new measurement values.
        std:: vector <std:: vector <float>> results;   // Results array accumulating observable values in time (will be printed to file).
        int k {0};                                     // Current index of results array to store measurements in.
        

        bool trajectory;    // If true, trajectory will be stored and printed to file.

        // Stores particle configurations in time (if --trajectory flag is set).
        std:: vector <std:: vector <coordinate<DIMENSION>>> trajectory_buffer;  
        
        std:: vector <float> times; // Times at which measurements are taken (printed to output file together with results).

        std:: vector <std:: string> col_names; // Names of the columns in the output file (names of the observables).

        const int N_iter;         // Number of interations in the simulation.
        const int N_meas;         // Take measurement any N_meas steps.
        
        const double stepsize;    // Stepsize used in the simulation. 
                                  // Needed here to obtain correct time values based on N_iter and N_meas for output file.         


        void add_to_results();

        // !!!!!!!!!!! DECLARE FUNCTIONS TO COMPUTE OBSERVABLES BASED ON SYSTEM'S CURRENT CONFIGURATION !!!!!!!!!!!!!!!!!!!!
        float get_center_of_mass_distance();
        float get_msd();
        float get_Tkin();

        // Add entry if desired.

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

};
// ##################### END OF CLASS DEFINITION ###########################################################################




// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################

template <size_t DIMENSION>
inline void measurement<DIMENSION>:: add_to_results(){
    /*
    Adds freshly computed observable values to results array storing
    time series that will be printed to file.
    */ 
    
    for (int i=0; i<no_observables; ++i) results[i][k] = observables[i];

    times[k] = k*N_meas*stepsize;   

    if (trajectory) trajectory_buffer[k] = model.positions;

    ++k;

}


template <size_t DIMENSION>
inline void measurement<DIMENSION>:: print_results(const std:: string outputname){
    /* 
    Print results to output.csv. Each observable will be printed to one column with 
    simulation time being the first column. 
    
    If bool trajectory is true, the whole trajectory will be printed to output.csv_trajectory 
    frame by frame, that is, the first N rows correspond to the first frame with the first column yielding
    the timestamp of that frame (identical for all N rows) and the next columns yielding a
    position coordinate each. The rows N+1 to 2N then give the next frame and so on.
    */

    // First write observable time series.
    std:: cout << "Writing results to file " << outputname << std:: endl;

    std:: ofstream file {outputname};  // Open file.
    
    // Write header with specified column names.
    file << "Time ";
    for (size_t k=0; k<col_names.size(); ++k){
        file << col_names[k] << " ";
    }
    file << "\n";

    // Write times and results array.
    for ( size_t i=0; i<times.size(); ++i )
    {
        file << times[i] << " "; 
        for (size_t j=0; j<no_observables; ++j) file << results[j][i] << " ";
        file << "\n";
    }

    file.close();


    // Write trajectory if needed.
    if (trajectory){
        std:: cout << "Writing trajectory...\n";
        std:: ofstream traj_file {outputname+"_trajectory"};
        
        // Write header with specified column names.
        traj_file << "Time";
        for (size_t dim=0; dim<DIMENSION; ++dim) traj_file << " dim" + std::to_string(dim+1);
        traj_file << "\n";

        // Write times and positions.
        double time;
        for ( size_t i=0; i<trajectory_buffer.size(); ++i )
        {
            time =  i*N_meas*stepsize;
            for ( size_t j=0; j<model.N_particles; ++j ){
                traj_file << time;
                for (size_t dim=0; dim<DIMENSION; ++dim){
                traj_file << " " << trajectory_buffer[i][j][dim];
                }
                traj_file << "\n";  
            }
        }

        traj_file.close();

    }
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPLEMENT OBSERVABLE CALCULATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

template <size_t DIMENSION>
inline float measurement<DIMENSION>:: get_center_of_mass_distance(){
/* 
    We use the method of  L. Bai and D. Breen, 
    ''Calculating Center of Mass in an Unbounded 2D Environment,'' 
     Journal of Graphics Tools, Vol. 13, No. 4, December 2008, pp. 53-60. 
*/

    // Compute center of mass (see the paper).
    const double two_L {2*model.L};
    const double pref {2*M_PI/model.L};
    const double pref2 {1/pref};
    coordinate<DIMENSION> center_of_mass, xi, zeta, theta;

    for (size_t n=0; n<model.N_particles; ++n){
        for (size_t dim=0; dim<DIMENSION; ++dim){
        
            theta[dim] = pref*model.positions[n][dim];
            xi[dim] += cos(theta[dim]);
            zeta[dim] += sin(theta[dim]);

        }
    }
    
    for (size_t dim=0; dim<DIMENSION; ++dim){
        xi[dim] *= pref2/model.N_particles;
        zeta[dim] *= pref2/model.N_particles;

        center_of_mass[dim] = pref2 * (atan2(-zeta[dim], -xi[dim]) + M_PI);
    }

    // Compute distance to COM.
    float dist {0};
    coordinate<DIMENSION> dist_dim;
    double sum;
    for (size_t n=0; n<model.N_particles; ++n){
        sum = 0;
        for (size_t dim=0; dim<DIMENSION; ++dim){
            dist_dim[dim] = model.positions[n][dim] - center_of_mass[dim];

            if (dist_dim[dim] > model.L)       dist_dim[dim] -= two_L;  // Periodic boundaries.
            else if (dist_dim[dim] < -model.L) dist_dim[dim] += two_L;
            
            sum += dist_dim[dim]*dist_dim[dim];
        }
        dist += sqrt(sum);
    }

    return dist/model.N_particles;

}    


template <size_t DIMENSION>
inline float measurement<DIMENSION>:: get_msd(){
    // Compute Mean Squared Displacement.

    coordinate<DIMENSION> diff;
    double msd {0}, two_L {2*model.L};

    for(int i=0; i<model.N_particles; ++i){
        
        for (size_t dim=0; dim<DIMENSION; ++dim){

            diff[dim] = model.positions[i][dim] - model.init_positions[i][dim];

            if (diff[dim] > model.L)         diff[dim] -= two_L;   // Periodic boundaries.
            else if (diff[dim] < -model.L)   diff[dim] += two_L;

            msd += diff[dim]*diff[dim];
        }
    }

    return 1./model.N_particles * msd;
}


template <size_t DIMENSION>
inline float measurement<DIMENSION>::  get_Tkin(){
    // Compute kinetic temperature.

    double Tkin {0};
    double v {0};

    for(int i=0; i<model.N_particles; ++i){
        for (size_t dim=0; dim<DIMENSION; ++dim){
            v = model.velocities[i][dim];
            Tkin += v*v;
        }
    }

    return Tkin/(2*model.N_particles);

}

// Add entries if needed.

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // MEASUREMENT_H