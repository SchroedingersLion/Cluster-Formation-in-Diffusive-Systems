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

constexpr double L {5};           // box volume = [-L,L]^n
int N_iter {12000};
constexpr int n_meas {1};
constexpr int n_part {1000};
constexpr double beta {300};  

constexpr double kappa {1./n_part};
constexpr int randomseed {1};
constexpr int THREADS = 10;

// ########## FORCE CONSTANTS ################

// Gaussian potential.
constexpr double sigma_2_gauss {1}; // sigma^2
constexpr double T_Tcrit_gauss = (1/beta) / (2*M_PI * n_part/(4*L*L) * sigma_2_gauss * 0.5*kappa);

// Morse potential.
constexpr double a_morse {1};
constexpr double r_morse {0.25};
constexpr double D_morse {4};
//############################################


constexpr double two_L {2*L};

// RNG.
std:: mt19937 twister;

struct coordinate{
    double x{0};
    double y{0};
};


static coordinate get_distances_ij(const coordinate position_i, const coordinate position_j){

    double dx {position_i.x - position_j.x};
    double dy {position_i.y - position_j.y};

    dx = dx > L ? dx - two_L : (dx < -L ? dx + two_L : dx);
    dy = dy > L ? dy - two_L : (dy < -L ? dy + two_L : dy);

    coordinate distances {dx, dy};

    return distances;

}

// ######## FORCES ###############################################################################

// Gaussian potential.
// static coordinate get_force_ij(const coordinate position_i, const coordinate position_j){

//     const coordinate dist {get_distances_ij(position_i, position_j)};  // gets (dx, dy) tupel.
    
//     const double dist_sq {dist.x*dist.x + dist.y*dist.y};

//     const double pref {-kappa/(2*sigma_2_gauss)};
//     const double expo_term {pref * exp(-dist_sq/(2*sigma_2_gauss))};

//     coordinate force_ij {expo_term*dist.x, expo_term*dist.y};

//     return force_ij;

// }


// Morse potential.
static coordinate get_force_ij(const coordinate position_i, const coordinate position_j){

    const coordinate dist {get_distances_ij(position_i, position_j)};  // gets (dx, dy) tupel.
    
    const double r {sqrt(dist.x*dist.x + dist.y*dist.y)};

    const double expo {exp(-a_morse * (r-r_morse))}; 
    const double pref {-kappa * a_morse * D_morse * (expo-expo*expo) / r};

    coordinate force_ij {pref*dist.x, pref*dist.y};

    return force_ij;

}


static std::vector<std::vector<coordinate>> forces_for_all_tasks( THREADS, std::vector<coordinate>(n_part, coordinate{ .x = 0.0, .y = 0.0 }));
static void compute_force_par(std::vector<coordinate>& forces, const std::vector<coordinate>& positions)
{

    for (auto& forces_for_specific_task : forces_for_all_tasks)
        std::fill(forces_for_specific_task.begin(), forces_for_specific_task.end(), coordinate{ .x = 0.0, .y = 0.0 });

    #pragma omp parallel for schedule(dynamic) num_threads(THREADS)
    for (int i = 0; i < n_part; ++i) {
        for (int j = i + 1; j < n_part; ++j) {

        coordinate force_ij = get_force_ij(positions[i], positions[j]);
        
        std::vector<coordinate>& forces_for_this_task = forces_for_all_tasks[omp_get_thread_num()];
        
        forces_for_this_task[i].x += force_ij.x;
        forces_for_this_task[i].y += force_ij.y;
        forces_for_this_task[j].x += -force_ij.x;
        forces_for_this_task[j].y += -force_ij.y;
        
        }
    }

    // Sum all of the task-specific forces into the output parameter.
    std::fill(forces.begin(), forces.end(), coordinate{ .x = 0.0, .y = 0.0 });
    for (auto const& forces_for_specific_task : forces_for_all_tasks)
        for (int i = 0; i < n_part; ++i)
        {
        forces[i].x += forces_for_specific_task[i].x;
        forces[i].y += forces_for_specific_task[i].y;
        }
}


static void U_step(std:: vector <coordinate>& positions, std:: vector <coordinate>& velocities, const double stepsize, const double gamma){

	std:: normal_distribution<> normal{0,1};

    // Prefactors used in the U step.
    const double pref_U1 {exp(-gamma*stepsize)};
    const double pref_U2 {(1-pref_U1)/gamma};
    const double pref_U3 {sqrt(2/(beta*gamma))};
    const double pref_U4 {sqrt(2*gamma/beta)};
    const double pref_Z1 {sqrt(stepsize)};
    const double pref_Z2 {sqrt( (1-pref_U1*pref_U1)/(2*gamma) )};
    const double pref_Z3 {sqrt( 2*(1-pref_U1)/(gamma*stepsize*(1+pref_U1)) )};
    const double pref_Z4 {sqrt( 1-2*(1-pref_U1)/(gamma*stepsize*(1+pref_U1)) )};

    const double pref_Z2_total1 {pref_Z2 * pref_Z3};
    const double pref_Z2_total2 {pref_Z2 * pref_Z4};


    // Compute new positions/velocities.
    coordinate Z1, Z2;
    coordinate xi1, xi2; 
    const size_t n_part {positions.size()};

    for (int i=0; i<n_part; ++i){

        xi1.x = normal(twister);
        xi1.y = normal(twister);
        xi2.x = normal(twister);
        xi2.y = normal(twister);

        Z1.x = pref_Z1 * xi1.x;
        Z1.y = pref_Z1 * xi1.y;

        Z2.x = pref_Z2_total1 * xi1.x + pref_Z2_total2 * xi2.x;
        Z2.y = pref_Z2_total1 * xi1.y + pref_Z2_total2 * xi2.y;

        positions[i].x += pref_U2 * velocities[i].x + pref_U3 * (Z1.x - Z2.x);
        positions[i].y += pref_U2 * velocities[i].y + pref_U3 * (Z1.y - Z2.y);

        velocities[i].x = pref_U1 * velocities[i].x + pref_U4 * Z2.x;
        velocities[i].y = pref_U1 * velocities[i].y + pref_U4 * Z2.y;

    }

    return;

}

static void B_step(std:: vector <coordinate>& velocities, std:: vector <coordinate>& forces, const double h){

    // Update velocities.
    for (int i=0; i<n_part; ++i){
        velocities[i].x += h*forces[i].x;
        velocities[i].y += h*forces[i].y;
    }

    return;

}

static void A_step(std:: vector <coordinate>& positions, const std:: vector <coordinate>& velocities, const double h){

        for (int i=0; i<n_part; ++i){
            positions[i].x += h*velocities[i].x;
            positions[i].y += h*velocities[i].y;
        }

    return;
}

static void O_step(std:: vector <coordinate>& velocities, const double h, const double gamma){

    std:: normal_distribution<> normal{0,1};
    const double a = exp(-gamma*h);
    const double pref = sqrt(1/beta *(1-a*a));

    for (int i=0; i<n_part; ++i){
        velocities[i].x = a*velocities[i].x + pref*normal(twister); 
        velocities[i].y = a*velocities[i].y + pref*normal(twister); 
    }

    return;
    
}


static void apply_periodic_boundaries(std:: vector <coordinate>& positions){

    double x, y;

    for (int i=0; i<n_part; ++i){
        x = positions[i].x;
        positions[i].x = x>L ? x-two_L : (x<-L ? x + two_L : x);
        y = positions[i].y;
        positions[i].y = y>L ? y-two_L : (y<-L ? y + two_L : y);
    }

    return;

}


static void print_to_csv(const std:: string& outputname, 
                        const std:: vector <std:: vector <coordinate>> positions,
                        const std:: vector <float> center_of_mass_distance,
                        const std:: vector <float> msds,
                        const std:: vector <float> Tkin,
                        const double stepsize, const double n_meas){    

    std:: ofstream file{outputname};
    double time;
    std:: cout << "Writing to file...\n";
    file << "Time  " << "x  " << "y  " << "com_distance " << "MSD " << "Tkin\n";
    for ( size_t i = 0; i<positions.size(); ++i )
    {
        time =  i*n_meas*stepsize;
        for ( size_t j = 0; j<n_part; ++j )
        {
            file << time << " " << positions[i][j].x << " " << positions[i][j].y << " " << center_of_mass_distance[i] << " " << msds[i] << " " << Tkin[i] << "\n";  
        }
    }

    file.close();

}


static void print_to_csv(const std:: string& outputname, 
                        const std:: vector <float> center_of_mass_distance,
                        const std:: vector <float> msds,
                        const std:: vector <float> Tkin,
                        const double stepsize, const double n_meas){    
// Overloaded version without printing positions.

    std:: ofstream file{outputname};
    double time;
    std:: cout << "Writing to file...\n";
    file << "Time " << "com_distance " << "MSD " << "Tkin\n";
    for ( size_t i = 0; i<center_of_mass_distance.size(); ++i )
    {
        time =  i*n_meas*stepsize;
        file << time << " " << center_of_mass_distance[i] << " " << msds[i] << " " << Tkin[i] << "\n";    
        
    }

    file.close();

}



static double get_center_of_mass_distance(const std:: vector <coordinate>& positions){
// We use the method of  
// L. Bai and D. Breen, ``Calculating Center of Mass in an Unbounded 2D Environment,'' 
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


static double get_msd(const std:: vector <coordinate>& positions, const std:: vector <coordinate>& init_positions){

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



static double get_Tkin(const std:: vector <coordinate>& velocities){

    double Tkin {0};

    for(auto vel : velocities){
        Tkin += vel.x*vel.x + vel.y*vel.y;
    }

    return Tkin/(2*velocities.size());

}


static std:: vector <coordinate> set_initial_positions(int seed){
    // uniform initialization.

    twister.seed(seed);
    std:: uniform_real_distribution<double> box_uniform(-L, L);
    std:: vector <coordinate> positions(n_part);
    for (int i=0; i<n_part; ++i){
        positions[i].x = box_uniform(twister);
        positions[i].y = box_uniform(twister);
    }
    
    return positions;

} 

static std:: vector <coordinate> set_initial_positions(){
    // square lattice initialization.

    std:: vector <coordinate> positions(n_part);
    
    // obtain number of particles per box dimension (for even spacing).
    int Nx = static_cast<int>(floor(sqrt(n_part)));
    while (n_part%Nx != 0) Nx += 1;
    int Ny = n_part/Nx;
    
    const double L_mod = L-0.00001; // to ensure no particle is placed on the edge.

    const double dx = 2*L_mod/Nx; // spacing between particles.
    const double dy = 2*L_mod/Ny;

    // place particles.
    for (int iy=0; iy<Ny; ++iy){
        for(int ix=0; ix<Nx; ++ix){

            positions[iy*Nx + ix].x = -L_mod + ix*dx;
            positions[iy*Nx + ix].y = -L_mod + iy*dy;

        }
    }

    return positions;

}


int main(int argc, char* argv[]){

    const double gamma = atof(argv[1]);
    const double h = atof(argv[2]);
    N_iter = atoi(argv[3]);
    int seed = atoi(argv[4]);

    // Output file name.
    std::ostringstream h_string;
    h_string << std::fixed << std::setprecision(2) << h;
    std:: string outputname {"results_h" + h_string.str() + "_gam" + argv[1] + "_N" + std::to_string(n_part) + "_seed" + std::to_string(seed) + ".csv"};
    
    std::cout<< "Simulation at T/Tcrit="<<T_Tcrit_gauss<<" with gamma="<<gamma<<std::endl;

    // Prepare simulation.

    // Set positions.
    // auto init_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const std:: vector <coordinate> init_positions = set_initial_positions(seed);
    // const std:: vector <coordinate> init_positions = set_initial_positions();
    std:: vector <coordinate> positions = init_positions;

    // Seed RNG for simulation.
    twister.seed(seed);

    // Set forces.
    std:: vector <coordinate> forces(n_part);
    compute_force_par(forces, positions);
    
    // Set velocities.
    std:: vector <coordinate> velocities(n_part);
    
    // Results vector.
    const size_t size_results_vector = N_iter / n_meas + 1;
    std:: vector <std:: vector <coordinate>> position_samples(size_results_vector);
    std:: vector <float>                     center_of_mass_distances(size_results_vector);
    std:: vector <float>                     msds(size_results_vector);
    std:: vector <float>                     Tkin(size_results_vector);
    int k {0};

    // Main loop.
    std:: normal_distribution<> normal{0,1};
    const double h_half = 0.5 * h;

    
    auto t1 = std:: chrono::high_resolution_clock::now();
    for (int i=0; i<=N_iter; ++i){
        
        if (i % n_meas == 0){  // Take measurement.
        
            position_samples[k] = positions;
            center_of_mass_distances[k] = get_center_of_mass_distance(positions);
            msds[k] = get_msd(positions, init_positions);
            Tkin[k] = get_Tkin(velocities);
            ++k;
        }

        // U_step(positions, velocities, h_half, gamma);
        // apply_periodic_boundaries(positions);
        // compute_force_par(forces, positions);
        // B_step(velocities, forces, h);
        // U_step(positions, velocities, h_half, gamma);

        B_step(velocities, forces, h_half);
        A_step(positions, velocities, h_half);
        apply_periodic_boundaries(positions);
        O_step(velocities, h, gamma);
        A_step(positions, velocities, h_half);
        apply_periodic_boundaries(positions);
        compute_force_par(forces, positions);
        B_step(velocities, forces, h_half);

        if (i%1000==0) std::cout << "Iteration "<<i<< " done!" << std::endl;

    }

    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";
    
    // print_to_csv(outputname, position_samples, center_of_mass_distances, msds, Tkin, h, n_meas);
    print_to_csv(outputname, center_of_mass_distances, msds, Tkin, h, n_meas);

    return 0;

}
