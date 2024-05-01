#include "IPS.h"

int main(int argc, char* argv[]){

    // Process command line arguments.
    auto [parse_result, options] = parseCommandLine(argc, argv); // Parse command line.
    if (parse_result.count("help")) {std:: cout << options.help() << std:: endl; return 0;}
    ParsedValues vals = processParsedValues(parse_result);


    // Set up simulation.
    IPS_model sys(vals.N_particles, vals.boxlength, vals.forcefield);

    measurement meas(vals.N_meas, vals.N_iter);

    simulation simu(sys, 
                    meas, 
                    vals.stepsize, 
                    vals.beta, 
                    vals.gamma, 
                    vals.N_iter, 
                    vals.threads, 
                    vals.integrator, 
                    vals.init_mode, 
                    vals.seed);

    // Run simulation.
    simu.run();

    // Print results.
    meas.print_to_csv(vals.output_name);


    return 0;

}