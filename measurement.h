#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <string>
#include <sstream>
#include <iomanip>

#include "IPS_model.h"

/*  This is the measurement class that is used by the simulation to take measurements of user-specified quantities. 
    In order to modify what observables are collected, the user has to do 2 things:
    a) in the constructor, they have to specify the numbers of observables to be taken by
       adjusting the value of the variable "no_observables". 
    b) in the function "take_measurement" the user has to adjust the formulas used to compute an observable from 
       the parameters. 
*/


// ###################### MEASUREMENT CLASS DEFINITION ##################################################

class measurement {

    public:

        // CONSTRUCTOR.
        measurement(const int method_type, const int burnin, const int t_meas, const int n_dist, const int max_iter)
            : n_meas {n_meas}, max_iter {max_iter}
            {
                
                /*######## ENTER THE NUMBER OF OBSERVABLES TO COLLECT ############*/
                no_observables = 3; 
                /*################################################################*/
                
                observables.resize(no_observables);
                results.resize(no_observables);
                col_names.resize(no_observables);
                
                int no_of_measurements {N_iter / n_meas + 1};
                
                for (auto& observable_vector : results) observable_vector.resize(no_of_measurements);
                times.resize(no_of_measurements);

            
            };



        void take_measurement(const IPS_model& model){

            /* ########### COMPUTE CURRENT OBSERVABLE VALUES FROM PARAMETERS ########
               The number of entries in vector "observables" must correspond to member variable "no_observables" set by the user
               in the constructor above. */          
            observables[0] = get_center_of_mass_distance(model);
            observables[1] = get_msd(model);
            observables[2] = get_Tkin(model);
            
            
            col_names[0] = "COM";
            col_names[1] = "MSD";
            col_names[2] = "Tkin";
            /*########################################################################*/

            add_to_results();  // Add new observables to results array and adds new time value.

            return;

        };


        void print_to_csv(const std:: string outputname);    // Prints out results array.


    private:
        int no_observables;                            // Number of observables to be taken.
        std:: vector <float> observables;              // Vector of size (no_observables) storing new measurement values.
        std:: vector <std:: vector <float>> results;   // Results array accumulating observable values in time (will be printed to file).
        int k {0};                                     // Current index of results array to store measurements in.
        
        float stepsize;            /* Stepsize used in the simulation (set by the simulation)
                                      This is needed to obtain the correct simulation times at which measurements are taken. */
        std:: vector <float> time; // Times at which measurements are taken (printed to output file together with results).

        std:: vector <std:: string> col_names; // Names of the columns in the output file (names of the observables).

        const int max_iter;  // Number of iterations in the simulation.
        const int n_meas;    // Number of iterations after which measure. (having these two in here is ugly actually).
        
        void add_to_results();

        float get_center_of_mass_distance(const IPS_model& model);
        float get_msd(const IPS_model& model);
        float get_Tkin(const IPS_model& model);

};
// ##################### END OF CLASS DEFINITION ##############################################




// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################


inline float measurement:: get_center_of_mass_distance(const IPS_model& model){
// We use the method of  L. Bai and D. Breen, 
// ''Calculating Center of Mass in an Unbounded 2D Environment,'' 
// Journal of Graphics Tools, Vol. 13, No. 4, December 2008, pp. 53-60. 

    // Compute center of mass.
    coordinate center_of_mass;
    const double pref = 2*M_PI/two_L;
    const double pref2 = 1/pref;
    coordinate xi, zeta, theta;

    for (const auto pos : positions){
        
        theta.x = pref*pos.x;
        xi.x += cos(theta.x);
        zeta.x += sin(theta.x);

        theta.y = pref*pos.y;
        xi.y += cos(theta.y);
        zeta.y += sin(theta.y);

    }
    xi.x *= pref2/n_part;
    xi.y *= pref2/n_part;
    zeta.x *= pref2/n_part;
    zeta.y *= pref2/n_part;

    center_of_mass.x = pref2 * (atan2(-zeta.x, -xi.x) + M_PI);
    center_of_mass.y = pref2 * (atan2(-zeta.y, -xi.y) + M_PI);


    // Compute distance to COM.
    float dist {0}, dist_x {0}, dist_y{0};
    for (const auto pos : positions){
        dist_x = pos.x - center_of_mass.x;
        dist_y = pos.y - center_of_mass.y;

        if (dist_x > L)       dist_x -= two_L;  // periodic boundaries.
        else if (dist_x < -L) dist_x += two_L;
        if (dist_y > L)       dist_y -= two_L;
        else if (dist_y < -L) dist_y += two_L;
        
        dist += sqrt(dist_x*dist_x + dist_y*dist_y);

    }

    return dist/n_part;

}    



inline float measurement:: get_msd(const IPS_model& model){
    
    coordinate diff;
    double msd {0};

    for(int i=0; i<positions.size(); ++i){
        
        diff.x = positions[i].x - init_positions[i].x;
        diff.y = positions[i].y - init_positions[i].y;

        if (diff.x > L)         diff.x -= two_L;
        else if (diff.x < -L)   diff.x += two_L;
        if (diff.y > L)         diff.y -= two_L;
        else if (diff.y < -L)   diff.y += two_L;

        msd += diff.x*diff.x + diff.y*diff.y;

    }

    return 1./positions.size() * msd;
}



inline float measurement:: get_Tkin(const IPS_model& model){

    double Tkin {0};

    for(auto vel : velocities){
        Tkin += vel.x*vel.x + vel.y*vel.y;
    }

    return Tkin/(2*velocities.size());

}



inline void measurement:: add_to_results(){

    for (int i=0; i<no_observables; ++i) results[i][k] = observables[i];
    times[k] = k*n_meas*stepsize;
    ++k;

}



inline void measurement:: print_to_csv(const std:: string outputname){

    // Overloaded version without printing positions.
    std:: cout << "Writing to file...\n";


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
        file << times[i] << " " 
        for (size_t j=0; j<no_observables; ++j) file << results[i][j] << " ";
        file << "\n";
    }

    file.close();

}


// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // MEASUREMENT_H