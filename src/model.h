#ifndef MODEL_H
#define MODEL_H

#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>


struct coordinate{    // Used to denote positions, velocities and forces.
    double x{0};
    double y{0};
};


// ###################### Model CLASS DEFINITION ##################################################################################

class IPS_model {

    public: 

        // MEMBERS THAT WILL BE ACCESSED BY MEASUREMENT OR SIMULATION CLASSES
        const double L;                // Box volume = [-L,L]^2.
        const int N_particles;         // No. of particles.
        std:: vector <coordinate> positions, init_positions, velocities, forces;
        coordinate (IPS_model::* get_force_ij) (const coordinate, const coordinate); // Points to interaction function between two particles.

        // CONSTRUCTOR.
        IPS_model(const int N_particles, const double boxlength, const std:: string& forcefield)
                : N_particles {N_particles}, 
                  L {boxlength/2}, 
                  kappa {1./N_particles}, 
                  forcefield {forcefield} 
            {

                std:: cout  << "Building IPS model with boxlength " << boxlength << " and " << N_particles << " particles.\n"
                            << "Forcefield: " << forcefield << std:: endl;

                // Resize vectors.
                positions.resize(N_particles);
                init_positions.resize(N_particles);
                velocities.resize(N_particles);
                forces.resize(N_particles);

                // Specify force field.
                if (forcefield=="gauss") get_force_ij = &IPS_model:: get_force_ij_gauss;
                else if (forcefield=="morse") get_force_ij = &IPS_model:: get_force_ij_morse;
                else throw std:: invalid_argument( "Invalid forcefield in model construction. Allowed are 'gauss' and 'morse'." );

            }

    private:

        const double kappa;             // Interaction potential prefactor (see paper).
        const std:: string forcefield;  // Specifies interaction potential.
        
        // ########## FORCES ###############################################################################################################
        coordinate get_distances_ij(const coordinate position_i, const coordinate position_j);
        
        // Gaussian potential.
        const double sigma_2_gauss {1}; // sigma^2.
        // const double T_Tcrit_gauss = (1/beta) / (2*M_PI * N_particles/(4*L*L) * sigma_2_gauss * 0.5*kappa);
        coordinate get_force_ij_gauss(const coordinate position_i, const coordinate position_j);            // Gauss interaction function.

        // Morse potential.
        const double a_morse {1};
        const double r_morse {0.25};
        const double D_morse {4};
        coordinate get_force_ij_morse(const coordinate position_i, const coordinate position_j);            // Morse interaction function.
        //##################################################################################################################################

};
// ##################### END OF CLASS DEFINITION ##############################################


// ##################### INLINE MEMBER FUNCTION DEFINITIONS ###################################
inline coordinate IPS_model:: get_distances_ij(const coordinate position_i, const coordinate position_j){

    const double two_L {2*L};

    double dx {position_i.x - position_j.x};
    double dy {position_i.y - position_j.y};

    dx = dx > L ? dx - two_L : (dx < -L ? dx + two_L : dx);
    dy = dy > L ? dy - two_L : (dy < -L ? dy + two_L : dy);

    coordinate distances {dx, dy};

    return distances;

}


// Gaussian potential.
inline coordinate IPS_model:: get_force_ij_gauss(const coordinate position_i, const coordinate position_j){

    const coordinate dist {get_distances_ij(position_i, position_j)};  // gets (dx, dy) tupel.
    
    const double dist_sq {dist.x*dist.x + dist.y*dist.y};

    const double pref {-kappa/(2*sigma_2_gauss)};
    const double expo_term {pref * exp(-dist_sq/(2*sigma_2_gauss))};

    coordinate force_ij {expo_term*dist.x, expo_term*dist.y};

    return force_ij;

}



// Morse potential.
inline coordinate IPS_model:: get_force_ij_morse(const coordinate position_i, const coordinate position_j){

    const coordinate dist {get_distances_ij(position_i, position_j)};  // Gets (dx, dy) tupel.
    
    const double r {sqrt(dist.x*dist.x + dist.y*dist.y)};

    const double expo {exp(-a_morse * (r-r_morse))}; 
    const double pref {-kappa * a_morse * D_morse * (expo-expo*expo) / r};

    coordinate force_ij {pref*dist.x, pref*dist.y};

    return force_ij;

}

// ########### END OF MEMBER DEFINITIONS ##############################################


#endif // MODEL