#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <chrono>
#include <omp.h>
#include <numbers>
#include <string>
#include <sstream>
#include <iomanip>

constexpr double L {5};           // box volume = [-L,L]^n
constexpr int N_iter {15000};
constexpr double h {0.1};
constexpr int n_meas {100};
constexpr int n_part {1000};
constexpr double beta {300};  

constexpr double kappa {1./n_part};
constexpr double sigma_2 {1};
constexpr int randomseed {2};

constexpr double h_half = 0.5 * h;
constexpr double T_Tcrit = (1/beta) / (2*std::numbers::pi * n_part/(4*L*L) * sigma_2 * 0.5*kappa);
constexpr double two_L {2*L};

// Prepare RNG.
std:: mt19937 twister;
std:: seed_seq seq{1,20,3200,403,5*randomseed+1,12000,73667,9474+randomseed,19151-randomseed};
std:: vector < std::uint32_t > seeds(1);


struct coordinate{
    double x{0};
    double y{0};
};


static double compute_sq_distances(const coordinate& pos1, const coordinate& pos2, 
                            double& dx, double& dy){
    
    dx = pos1.x - pos2.x;
    if (dx > L)       dx -= two_L;
    else if (dx < -L) dx += two_L;


    dy = pos1.y - pos2.y;
    if (dy > L)       dy -= two_L;
    else if (dy < -L) dy += two_L;
    
    return dx*dx + dy*dy;

}

static double force_term(double& dist2){

     double pref {-kappa/(2*sigma_2)};
     double expo;
     expo = exp(-dist2/(2*sigma_2));
     
    return pref*expo;

}

static void compute_force(std:: vector<coordinate>& forces, const std:: vector<coordinate>& positions){

     double dx, dy;
     double dist2;
     double force;

    for (int i=0; i<n_part; ++i){
        forces[i].x = forces[i].y = 0;
    }

    for (int i=0; i<n_part; ++i){
        for(int j=i+1; j<n_part; ++j){
            dist2 = compute_sq_distances(positions[i], positions[j], dx, dy);
            force = force_term(dist2);
            forces[i].x += force*dx;
            forces[i].y += force*dy;
            forces[j].x += -force*dx;
            forces[j].y += -force*dy;
        }
    }

    return;

}

constexpr int THREADS = 10;
static std::vector<std::vector<coordinate>> forces_for_all_tasks( THREADS, std::vector<coordinate>(n_part, coordinate{ .x = 0.0, .y = 0.0 }));
static void compute_force_par(std::vector<coordinate>& forces, const std::vector<coordinate>& positions)
{


double pref {-kappa/(2*sigma_2)};
for (auto& forces_for_specific_task : forces_for_all_tasks)
    std::fill(forces_for_specific_task.begin(), forces_for_specific_task.end(), coordinate{ .x = 0.0, .y = 0.0 });

#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
  for (int i = 0; i < n_part; ++i) {
    for (int j = i + 1; j < n_part; ++j) {

      double dx = positions[i].x - positions[j].x;
      dx = dx > L ? dx - two_L : (dx < -L ? dx + two_L : dx);

      double dy = positions[i].y - positions[j].y;
      dy = dy > L ? dy - two_L : (dy < -L ? dy + two_L : dy);

      double dist2 = dx * dx + dy * dy;

      double expo = exp(-dist2 / (2 * sigma_2));

      double force = pref * expo;

      std::vector<coordinate>& forces_for_this_task = forces_for_all_tasks[omp_get_thread_num()];
      forces_for_this_task[i].x += force * dx;
      forces_for_this_task[i].y += force * dy;
      forces_for_this_task[j].x += -force * dx;
      forces_for_this_task[j].y += -force * dy;
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
    const double pref_U3 {sqrt(2/beta*gamma)};
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
                  const double stepsize, const double n_meas){    

    std:: ofstream file{outputname};
    double time;
    std:: cout << "Writing to file...\n";
    file << "Time  " << "x  " << "y  " << "com_distance\n";
    for ( size_t i = 0; i<positions.size(); ++i )
    {
        time =  i*n_meas*stepsize;
        for ( size_t j = 0; j<positions[0].size(); ++j )
        {
            file << time << " " << positions[i][j].x << " " << positions[i][j].y << " " << center_of_mass_distance[i] << "\n";  
        }
    }

    file.close();

}


static double get_center_of_mass_distance(const std:: vector <coordinate>& positions){

    coordinate center_of_mass;
    for (const auto pos : positions){
        center_of_mass.x += pos.x;
        center_of_mass.y += pos.y;
    }

    center_of_mass.x /= n_part;
    center_of_mass.y /= n_part;

    float dist {0}, dist_x {0}, dist_y{0};
    for (const auto pos : positions){
        dist_x = pos.x - center_of_mass.x;
        dist_y = pos.y - center_of_mass.y;

        if (dist_x > L)       dist_x -= two_L;  // periodic boundaries
        else if (dist_x < -L) dist_x += two_L;
        if (dist_y > L)       dist_y -= two_L;
        else if (dist_y < -L) dist_y += two_L;

        dist += sqrt(dist_x*dist_x + dist_y*dist_y);
    }

    return dist/n_part;

}



int main(int argc, char* argv[]){

    double gamma = atof(argv[1]);
    std::ostringstream h_string;
    h_string << std::fixed << std::setprecision(1) << h;
    std:: string outputname {"trajectory_h" + h_string.str() + "_gam" + argv[1] +".csv"};

    std::cout<< "Simulation at T/Tcrit="<<T_Tcrit<<" with gamma="<<gamma<<std::endl;

    // Prepare simulation.
    seq.generate(seeds.begin(), seeds.end());
    twister.seed(seeds.at(0)); 

    // Set positions.
    std:: uniform_real_distribution<double> box_uniform(-L, L);
    std:: vector <coordinate> positions(n_part);
    for (int i=0; i<n_part; ++i){
        positions[i].x = box_uniform(twister);
        positions[i].y = box_uniform(twister);
    }

    // Set forces.
    std:: vector <coordinate> forces(n_part);
    compute_force_par(forces, positions);
    
    // Set velocities.
    std:: vector <coordinate> velocities(n_part);
    
    // Results vector.
    std:: vector < std:: vector <coordinate>> position_samples(N_iter / n_meas + 1);
    std:: vector <float> center_of_mass_distances(N_iter / n_meas + 1);
    int k {0};

    // Main loop.
    std:: normal_distribution<> normal{0,1};
    
    auto t1 = std:: chrono::high_resolution_clock::now();
    for (int i=0; i<=N_iter; ++i){
        
        if (i % n_meas == 0){  // Take measurement.
        
            position_samples[k] = positions;
            center_of_mass_distances[k] = get_center_of_mass_distance(positions);
            ++k;
        }

        U_step(positions, velocities, h_half, gamma);
        apply_periodic_boundaries(positions);
        compute_force_par(forces, positions);
        B_step(velocities, forces, h);
        U_step(positions, velocities, h_half, gamma);

        // B_step(velocities, forces, h_half);
        // A_step(positions, velocities, h_half);
        // apply_periodic_boundaries(positions);
        // O_step(velocities, h, gamma);
        // A_step(positions, velocities, h_half);
        // apply_periodic_boundaries(positions);
        // compute_force_par(forces, positions);
        // B_step(velocities, forces, h_half);

        if (i%1000==0) std::cout << "Iteration "<<i<< " done!\n";

    }
    
    auto t2 = std:: chrono:: high_resolution_clock:: now();
    auto ms_int = std:: chrono:: duration_cast < std:: chrono:: seconds > (t2 - t1);
    std:: cout << "Execution took " << ms_int.count() << " seconds!\n";
    
    print_to_csv(outputname, position_samples, center_of_mass_distances, h, n_meas);

    return 0;

}