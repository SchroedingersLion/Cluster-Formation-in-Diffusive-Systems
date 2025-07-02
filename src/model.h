#ifndef MODEL_H
#define MODEL_H

#include "coordinate.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>

#define _USE_MATH_DEFINES

/*  
Holds the model class which specifies a cubic simulation box with N particles that interact via a pairwise
interaction potential. In particular, this class implements the pairwise forces given by various interaction potentials.
It is straightforward to extend the code by adding new pairwise interactions.

All member functions are implemented inline for simplicity.
*/



// ###################### Model CLASS DEFINITION ##################################################################################
template <size_t DIMENSION>
class IPS_model {

    public: 

        // MEMBERS THAT WILL BE ACCESSED BY MEASUREMENT OR SIMULATION CLASSES.
        const double L;                // Box volume = [-L,L]^2.
        const int N_particles;         // No. of particles.
        
        std:: vector <coordinate<DIMENSION>> positions, init_positions, velocities, forces;
        
        // Point to interaction function between two particles.
        void (IPS_model::* get_force_ij) (const coordinate<DIMENSION>&, coordinate<DIMENSION>&); 

        // Get (dx, dy) tupel of distances between particles i and j.
        void get_distances_ij(const size_t i, const size_t j, coordinate<DIMENSION>& distances);

        // DEFAULT CONSTRUCTOR.
        IPS_model<DIMENSION>(){};

        // CONSTRUCTOR.
        IPS_model<DIMENSION>(const int N_particles, const double boxlength, const std:: string& forcefield)
                : N_particles {N_particles}, 
                  L {boxlength/2}, 
                  kappa {1./N_particles}, 
                  forcefield {forcefield}
            {

                std:: cout  << "Building IPS model with boxlength " << boxlength << " and " << N_particles << " particles.\n"
                            << "Forcefield: " << forcefield << std:: endl;

                // Resize vectors.
                resize_vectors();

                // Specify force field.
                if (forcefield=="gauss") get_force_ij = &IPS_model:: get_force_ij_gauss;
                else if (forcefield=="morse") get_force_ij = &IPS_model:: get_force_ij_morse;
                else if (forcefield=="gem4") get_force_ij = &IPS_model:: get_force_ij_gem4;
                else throw std:: invalid_argument( "Invalid forcefield in model construction. Allowed are 'gauss', 'morse', and 'gem4'." );

            }

    private:

        const double kappa;             // Interaction potential prefactor (see paper).
        const std:: string forcefield;  // Specifies interaction potential.
        void resize_vectors();
        
        // ########## FORCES ###############################################################################################################        

        void get_force_ij_gauss(const coordinate<DIMENSION>& distance, coordinate<DIMENSION>& force_ij);  // Gauss interaction function.

        void get_force_ij_morse(const coordinate<DIMENSION>& distance, coordinate<DIMENSION>& force_ij);  // Morse interaction function.

        void get_force_ij_gem4(const coordinate<DIMENSION>& distance, coordinate<DIMENSION>& force_ij);   // GEM-4 interaction function.
        
        double my_pow(double x, int n); // Help function to compute x^n.

        //##################################################################################################################################

};
// ##################### END OF CLASS DEFINITION ##############################################



// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################
template <size_t DIMENSION>
inline void IPS_model<DIMENSION>:: resize_vectors(){

    positions.resize(N_particles, coordinate<DIMENSION>());
    init_positions.resize(N_particles, coordinate<DIMENSION>());
    velocities.resize(N_particles, coordinate<DIMENSION>());
    forces.resize(N_particles, coordinate<DIMENSION>());

}


template <size_t DIMENSION>
inline void IPS_model<DIMENSION>:: get_distances_ij(const size_t i, const size_t j, coordinate<DIMENSION>& distances){

    const double two_L {2*L};
    double dx;

    for (size_t dim = 0; dim<DIMENSION; ++dim){
        dx = positions[i][dim] - positions[j][dim];
        dx = dx > L ? dx - two_L : (dx < -L ? dx + two_L : dx);  // Correct for periodic boundaries.
        distances[dim] = dx;
    }

    return;

}


// Gaussian interaction.
const double sigma_2_gauss {0.5}; // sigma^2.
template <size_t DIMENSION>
inline void IPS_model<DIMENSION>:: get_force_ij_gauss(const coordinate<DIMENSION>& distance, coordinate<DIMENSION>& force_ij){

    
    double dist_sq {0};
    for (size_t dim = 0; dim<DIMENSION; ++dim) dist_sq += distance[dim]*distance[dim];

    const double pref {-kappa/(sigma_2_gauss)};
    const double expo_term {pref * exp(-dist_sq/(2*sigma_2_gauss))};

    for (size_t dim = 0; dim<DIMENSION; ++dim) force_ij[dim] = expo_term*distance[dim];

    return;

}


// GEM-Alpha interaction (reduces to Gaussian for alpha=2, but is less efficient compute-wise).
const double sigma_2_gem4 {0.5};
const double sqrt_two_sigma_2_gem4 {sqrt(2*sigma_2_gem4)};
const int alpha_gem4 {4};
const double pref_gem4 {alpha_gem4 / sqrt_two_sigma_2_gem4};

template <size_t DIMENSION>
inline void IPS_model<DIMENSION>:: get_force_ij_gem4(const coordinate<DIMENSION>& distance, coordinate<DIMENSION>& force_ij){

    double dist {0};
    for (size_t dim = 0; dim<DIMENSION; ++dim) dist += distance[dim]*distance[dim];
    dist = sqrt(dist);

    const double normed_distance {dist/ sqrt_two_sigma_2_gem4};
    const double exponent {my_pow(normed_distance, alpha_gem4)};

    const double pref {-kappa*pref_gem4*exponent/normed_distance};
    const double expo_term {pref * exp(-exponent)};
    
    for (size_t dim = 0; dim<DIMENSION; ++dim) force_ij[dim] =  distance[dim] > 0 ? expo_term : -expo_term;
    

    return;

}


template <size_t DIMENSION>
inline double IPS_model<DIMENSION>:: my_pow(double x, int n){
    double r = 1.0;

    while(n > 0){
        r *= x;
        --n;
    }

    return r;
}



// Morse potential.
const double a_morse {2};
const double r_morse {0}; 
const double D_morse {1};

template <size_t DIMENSION>
inline void IPS_model<DIMENSION>:: get_force_ij_morse(const coordinate<DIMENSION>& distance, coordinate<DIMENSION>& force_ij){
    
    double dist {0};
    for (size_t dim = 0; dim<DIMENSION; ++dim) dist += distance[dim]*distance[dim];
    dist = sqrt(dist);

    const double expo {exp(-a_morse * (dist-r_morse))}; 
    const double pref {kappa * 2* a_morse * D_morse * (expo*expo-expo) / dist};

    for (size_t dim = 0; dim<DIMENSION; ++dim) force_ij[dim] =  pref*distance[dim];

    return;

}

// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // MODEL_H