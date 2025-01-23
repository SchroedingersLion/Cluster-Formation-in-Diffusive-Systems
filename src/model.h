#ifndef MODEL_H
#define MODEL_H

#include "coordinate.h"

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>




// ###################### Model CLASS DEFINITION ##################################################################################

class IPS_model {

    public: 

        // MEMBERS THAT WILL BE ACCESSED BY MEASUREMENT OR SIMULATION CLASSES
        const double L;                // Box volume = [-L,L]^2.
        const int N_particles;         // No. of particles.
        const size_t dimension;
        std:: vector <coordinate> positions, init_positions, velocities, forces;
        void (IPS_model::* get_force_ij) (const coordinate&, coordinate&); // Points to interaction function between two particles.
        void get_distances_ij(const size_t i, const size_t j, coordinate& distances);                 // Get (dx, dy) tupel of distances between particles i,j.

        // CONSTRUCTOR.
        IPS_model(const int N_particles, const double boxlength, const std:: string& forcefield, const size_t dimension)
                : N_particles {N_particles}, 
                  L {boxlength/2}, 
                  kappa {1./N_particles}, 
                  forcefield {forcefield},
                  dimension {dimension} 
            {

                std:: cout  << "Building IPS model with boxlength " << boxlength << " and " << N_particles << " particles.\n"
                            << "Forcefield: " << forcefield << std:: endl;

                // Resize vectors.
                resize_vectors();

                // Specify force field.
                if (forcefield=="gauss") get_force_ij = &IPS_model:: get_force_ij_gauss;
                // else if (forcefield=="morse") get_force_ij = &IPS_model:: get_force_ij_morse;
                // else if (forcefield=="gem4") get_force_ij = &IPS_model:: get_force_ij_gem4;
                else throw std:: invalid_argument( "Invalid forcefield in model construction. Allowed are 'gauss', 'morse', and 'gem4'." );

            }

    private:

        const double kappa;             // Interaction potential prefactor (see paper).
        const std:: string forcefield;  // Specifies interaction potential.
        void resize_vectors();
        
        // ########## FORCES ###############################################################################################################        
        // Gaussian potential.
        void get_force_ij_gauss(const coordinate& distance, coordinate& force_ij);            // Gauss interaction function.

        // Morse potential.
        // coordinate get_force_ij_morse(const coordinate distance);            // Morse interaction function.

        // // GEM-4 potential
        // coordinate get_force_ij_gem4(const coordinate distance);            // GEM-4 interaction function.
        
        double my_pow(double x, int n); // Help function to compute x^n.

        //##################################################################################################################################

};
// ##################### END OF CLASS DEFINITION ##############################################


// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################
inline void IPS_model:: resize_vectors(){

    positions.resize(N_particles, coordinate(dimension));
    init_positions.resize(N_particles, coordinate(dimension));
    velocities.resize(N_particles, coordinate(dimension));
    forces.resize(N_particles, coordinate(dimension));

}


inline void IPS_model:: get_distances_ij(const size_t i, const size_t j, coordinate& distances){

    const double two_L {2*L};
    double dx;

    for (size_t dim = 0; dim<dimension; ++dim){
        dx = positions[i][dim] - positions[j][dim];
        dx = dx > L ? dx - two_L : (dx < -L ? dx + two_L : dx);
        distances[dim] = dx;
    }

    return;

}


// Gaussian potential.
const double sigma_2_gauss {0.5}; // sigma^2.
inline void IPS_model:: get_force_ij_gauss(const coordinate& distance, coordinate& force_ij){

    
    double dist_sq {0};
    for (size_t dim = 0; dim<dimension; ++dim) dist_sq += distance[dim]*distance[dim];

    const double pref {-kappa/(sigma_2_gauss)};
    const double expo_term {pref * exp(-dist_sq/(2*sigma_2_gauss))};

    for (size_t dim = 0; dim<dimension; ++dim) force_ij[dim] = expo_term*distance[dim];

    return;

}

// // Gaussian potential.
// const double sigma_gauss {sqrt(0.5)}; // Sigma.
// const int alpha_gauss {2};         // Exponent in exponential (==2 for Gauss).

// const double prefactor_gauss {alpha_gauss/pow(sqrt(2)*sigma_gauss, alpha_gauss)};   // Help constants.
// const double denominator_gauss {sqrt(2)*sigma_gauss};

// inline coordinate IPS_model:: get_force_ij_gauss(const coordinate position_i, const coordinate position_j){

//     const coordinate dist {get_distances_ij(position_i, position_j)};  // Get (dx, dy) tupel.
    
//     const double distance {sqrt(dist.x*dist.x + dist.y*dist.y)};

//     const double expo_term {-kappa * prefactor_gauss * my_pow(distance, alpha_gauss-2) * exp( - my_pow(distance/denominator_gauss, alpha_gauss) )};
//     // const double expo_term {-kappa * prefactor_gauss  * exp( - my_pow(distance/denominator_gauss, alpha_gauss) )};


//     coordinate force_ij {expo_term*dist.x, expo_term*dist.y};

//     return force_ij;

// }


// GEM-Alpha potential (reduces to Gaussian for alpha=2, but is less efficient compute-wise).
// const double sigma_2_gem4 {0.5};
// const double sqrt_two_sigma_2_gem4 {sqrt(2*sigma_2_gem4)};
// const int alpha_gem4 {4};
// const double pref_gem4 {alpha_gem4 / sqrt_two_sigma_2_gem4};
// inline coordinate IPS_model:: get_force_ij_gem4(const coordinate distance){

//     const double normed_distance {sqrt(distance.x*distance.x) / sqrt_two_sigma_2_gem4};
//     const double exponent {my_pow(normed_distance, alpha_gem4)};

//     const double pref {-kappa*pref_gem4*exponent/normed_distance};
//     const double expo_term {pref * exp(-exponent)};

//     coordinate force_ij {distance.x > 0 ? expo_term : -expo_term};

//     return force_ij;

// }


// inline double IPS_model:: my_pow(double x, int n){
//     double r = 1.0;

//     while(n > 0){
//         r *= x;
//         --n;
//     }

//     return r;
// }


// // Morse potential.
// const double a_morse {2};
// const double r_morse {0}; 
// const double D_morse {1};

// inline coordinate IPS_model:: get_force_ij_morse(const coordinate distance){
    
//     const double dist {sqrt(distance.x*distance.x)};

//     const double expo {exp(-a_morse * (dist-r_morse))}; 
//     const double pref {kappa * 2* a_morse * D_morse * (expo*expo-expo) / dist};

//     coordinate force_ij {pref*distance.x};

//     return force_ij;

// }

// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // MODEL_H