#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>

#include "model.h"



// ###################### MEASUREMENT CLASS DEFINITION ##################################################

class measurement {

    public:

        // CONSTRUCTOR.
        measurement(const IPS_model& model, int N_meas, const int N_iter, const double stepsize, const bool trajectory)
            : model {model}, N_meas {N_meas}, N_iter {N_iter}, stepsize {stepsize}, trajectory {trajectory}
            {
                
                /*######## ENTER THE NUMBER OF OBSERVABLES TO COLLECT ############*/
                no_observables = 3; 
                /*################################################################*/
                
                observables.resize(no_observables);
                results.resize(no_observables);
                col_names.resize(no_observables);

                /*################# ENTER NAMES OF OBSERVABLES (WILL BE HEADER OF OUTPUTFILE)####*/
                col_names[0] = "COM";
                col_names[1] = "MSD";
                col_names[2] = "Tkin";
                /*################################################################################*/
                
                int no_of_measurements {N_iter / N_meas + 1};
                
                for (auto& observable_vector : results) observable_vector.resize(no_of_measurements);
                times.resize(no_of_measurements);
                if (trajectory) trajectory_buffer.resize(no_of_measurements);

            
            };



        void take_measurement(){

            /* ########### COMPUTE CURRENT OBSERVABLE VALUES FROM PARAMETERS ########
               The number of entries in vector "observables" must correspond to member variable "no_observables" set by the user
               in the constructor above. */          
            
            observables[0] = get_center_of_mass_distance();
            observables[1] = get_msd();
            observables[2] = get_Tkin();
            /*########################################################################*/

            add_to_results();  // Add new observables to results array and add new time value.

            return;

        };


        void print_results(const std:: string outputname);    // Prints results array to .csv.


    private:
        const IPS_model& model;                              // Model to take measurements on.
        int no_observables;                            // Number of observables to be taken.
        std:: vector <float> observables;              // Vector of size (no_observables) storing new measurement values.
        std:: vector <std:: vector <float>> results;   // Results array accumulating observable values in time (will be printed to file).
        int k {0};                                     // Current index of results array to store measurements in.
        
        bool trajectory;    // If true, trajectory will stored and printed to file.
        std:: vector <std:: vector <coordinate>> trajectory_buffer;  // Stores particle configurations in time (if --trajectory flag is set).
        
        std:: vector <float> times; // Times at which measurements are taken (printed to output file together with results).

        std:: vector <std:: string> col_names; // Names of the columns in the output file (names of the observables).

        const int N_iter;  // Number of interations in the simulation.
        const int N_meas;  // Take measurement any N_meas steps.
        const double stepsize;    // Stepsize used in the simulation        
                                      
        void add_to_results();

        float get_center_of_mass_distance();
        float get_msd();
        float get_Tkin();

};
// ##################### END OF CLASS DEFINITION ##############################################




// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################


inline float measurement:: get_center_of_mass_distance(){
// We use the method of  L. Bai and D. Breen, 
// ''Calculating Center of Mass in an Unbounded 2D Environment,'' 
// Journal of Graphics Tools, Vol. 13, No. 4, December 2008, pp. 53-60. 

    // Compute center of mass.
    coordinate center_of_mass;
    const double two_L {2*model.L};
    const double pref {2*M_PI/two_L};
    const double pref2 {1/pref};
    coordinate xi, zeta, theta;

    for (const auto pos : model.positions){
        
        theta.x = pref*pos.x;
        xi.x += cos(theta.x);
        zeta.x += sin(theta.x);

        theta.y = pref*pos.y;
        xi.y += cos(theta.y);
        zeta.y += sin(theta.y);

    }
    xi.x *= pref2/model.N_particles;
    xi.y *= pref2/model.N_particles;
    zeta.x *= pref2/model.N_particles;
    zeta.y *= pref2/model.N_particles;

    center_of_mass.x = pref2 * (atan2(-zeta.x, -xi.x) + M_PI);
    center_of_mass.y = pref2 * (atan2(-zeta.y, -xi.y) + M_PI);


    // Compute distance to COM.
    float dist {0}, dist_x {0}, dist_y{0};
    for (const auto pos : model.positions){
        dist_x = pos.x - center_of_mass.x;
        dist_y = pos.y - center_of_mass.y;

        if (dist_x > model.L)       dist_x -= two_L;  // periodic boundaries.
        else if (dist_x < -model.L) dist_x += two_L;
        if (dist_y > model.L)       dist_y -= two_L;
        else if (dist_y < -model.L) dist_y += two_L;
        
        dist += sqrt(dist_x*dist_x + dist_y*dist_y);

    }

    return dist/model.N_particles;

}    



inline float measurement:: get_msd(){
    
    coordinate diff;
    double msd {0}, two_L {2*model.L};

    for(int i=0; i<model.positions.size(); ++i){
        
        diff.x = model.positions[i].x - model.init_positions[i].x;
        diff.y = model.positions[i].y - model.init_positions[i].y;

        if (diff.x > model.L)         diff.x -= two_L;
        else if (diff.x < -model.L)   diff.x += two_L;
        if (diff.y > model.L)         diff.y -= two_L;
        else if (diff.y < -model.L)   diff.y += two_L;

        msd += diff.x*diff.x + diff.y*diff.y;

    }

    return 1./model.positions.size() * msd;
}



inline float measurement:: get_Tkin(){

    double Tkin {0};

    for(auto vel : model.velocities){
        Tkin += vel.x*vel.x + vel.y*vel.y;
    }

    return Tkin/(2*model.velocities.size());

}



inline void measurement:: add_to_results(){

    for (int i=0; i<no_observables; ++i) results[i][k] = observables[i];
    times[k] = k*N_meas*stepsize;
    if (trajectory) trajectory_buffer[k] = model.positions;
    ++k;

}



inline void measurement:: print_results(const std:: string outputname){

    std:: cout << "Writing results to file..." << std:: endl;

    std:: ofstream file {outputname};
    
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
        traj_file << "Time " << "x " << "y\n";

        // Write times and positions.
        size_t number_particles {trajectory_buffer[0].size()};
        double time;
        for ( size_t i=0; i<trajectory_buffer.size(); ++i )
        {
            time =  i*N_meas*stepsize;
            for ( size_t j=0; j<number_particles; ++j )
            {
                traj_file << time << " " << trajectory_buffer[i][j].x << " " << trajectory_buffer[i][j].y <<  "\n";  
            }
        }

        traj_file.close();
    }


}



// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // MEASUREMENT_H