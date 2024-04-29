#include "IPS.h"

int main(int argc, char* argv[]){

    const double gamma = atof(argv[1]);
    const double h = atof(argv[2]);
    N_iter = atoi(argv[3]);
    int seed = atoi(argv[4]);

    int n_meas = 1;

    double stepsize=0.01; 
    double beta=300; 
    double gamma=1; 
    int Niter=10000;

    std:: string outputname = "test_output.csv";

    IPS_model sys();

    measurement meas(n_meas, N_iter);

    simulation simu(sys, meas, stepsize, beta, gamma, Niter);

    simu.run();

    meas.print_to_csv(outputname);

    return 0;

}