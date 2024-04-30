#include "IPS.h"

int main(int argc, char* argv[]){


    // Parse command line
    auto [result, options] = parseCommandLine(argc, argv);

    // Check if help option is provided
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    // Process parsed values
    ParsedValues values = processParsedValues(result);




    // const double gamma = atof(argv[1]);
    // const double h = atof(argv[2]);
    // const int N_iter = atoi(argv[3]);
    // const int seed = atoi(argv[4]);

    // int n_meas = 1;

    // double beta=300; 

    // std:: string outputname = "test_output.csv";

    // std::string force = "morse";
    
    IPS_model sys(1000, 5, force);

    measurement meas(n_meas, N_iter);

    simulation simu(sys, meas, h, beta, gamma, N_iter);

    simu.run();

    meas.print_to_csv(outputname);



    return 0;

}