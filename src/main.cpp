#include "IPS.h"


int main(int argc, char* argv[]){

    // Process command line arguments.
    auto [parse_result, options] = parseCommandLine(argc, argv); 
    if (parse_result.count("help")) {std:: cout << options.help() << std:: endl; return 0;}
    ParsedValues vals = processParsedValues(parse_result);


    // Set up simulation. Switch statement for dimensionality.
    switch(vals.dimension){
        
        case 1: {
            IPS_model<1> sys(vals.N_particles, vals.boxlength, vals.forcefield);
            simulation<1> simu(sys, 
                    vals.stepsize, 
                    vals.beta, 
                    vals.gamma, 
                    vals.N_iter,
                    vals.N_meas, 
                    vals.threads, 
                    vals.integrator, 
                    vals.init_mode, 
                    vals.seed,
                    vals.trajectory);
            
            simu.run();
            simu.meas.print_results(vals.output_name);

            break;
        }
        
        case 2: {
            IPS_model<2> sys(vals.N_particles, vals.boxlength, vals.forcefield);
            simulation<2> simu(sys, 
                    vals.stepsize, 
                    vals.beta, 
                    vals.gamma, 
                    vals.N_iter,
                    vals.N_meas, 
                    vals.threads, 
                    vals.integrator, 
                    vals.init_mode, 
                    vals.seed,
                    vals.trajectory);
            
            simu.run();
            simu.meas.print_results(vals.output_name);

            break;
        }

        case 3: {
            IPS_model<3> sys(vals.N_particles, vals.boxlength, vals.forcefield);
            simulation<3> simu(sys, 
                    vals.stepsize, 
                    vals.beta, 
                    vals.gamma, 
                    vals.N_iter,
                    vals.N_meas, 
                    vals.threads, 
                    vals.integrator, 
                    vals.init_mode, 
                    vals.seed,
                    vals.trajectory);
            
            simu.run();
            simu.meas.print_results(vals.output_name);

            break;
        }

        default:
            throw std::out_of_range("\nPassed inadmissible dimension argument. Allowed is `1`, `2`, or `3`.\n");
            break;
    
    }
    

    
    return 0;

}
