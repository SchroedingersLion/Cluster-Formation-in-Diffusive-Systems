#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <string>
#include <iostream>
#include <cxxopts.hpp>

    gamma
    stepsize
    N_iter
    beta
    seed
    meas_freq
    outputfile_name
    forcefield
    N_part
    threads
    integrator
    init_mode
    boxlength



std::pair<cxxopts::ParseResult, cxxopts::Options> parseCommandLine(int argc, char* argv[]) {
    cxxopts::Options options("MyProgram", "Description");

    // Define command line options
    options.add_options()
        ("Nparticles",  "Number of particles.",                                         cxxopts::value<int>()->default_value("1000"))
        ("boxlength",   "Length of edge of square simulation box.",                     cxxopts::value<double>()->default_value("10"))
        ("forcefield",  "Forcefield between two particles. Either 'gauss' or 'morse'.", cxxopts::value<std::string>()->default_value("gauss"))
        ("init_mode",   "Initial positions of the system. Either 'uniform' or 'grid'.", cxxopts::value<std::string>()->default_value("uniform"))
        ("beta",        "Inverse temperature parameter in Langevin dynamics.",          cxxopts::value<double>()->default_value("10"))
        ("gamma",       "Friction parameter in Langevin dynamics.",                     cxxopts::value<double>()->default_value("0.1"))
        ("Niter",       "Number of simulation steps.",                                  cxxopts::value<int>()->default_value("10000"))
        ("stepsize",    "Simulation stepsize.",                                         cxxopts::value<double>()->default_value("0.01"))
        ("Nmeas",       "Take a measurement any 'Nmeas' iterations.",                   cxxopts::value<int>()->default_value("10"))
        ("integrator",  "Integrator to be used. Either 'BAOAB' or 'UBU'.",              cxxopts::value<std::string>()->default_value("BAOAB"))
        ("threads",     "Number of threads used in the force evaluation.",              cxxopts::value<int>()->default_value("4"))
        ("seed",        "Randomseed.",                                                  cxxopts::value<int>()->default_value("1"))
        ("output_name", "Name of the printed file holding the results.",                cxxopts::value<std::string>()->default_value("results_test.csv"))               
        ("help",        "Print help");

    // Parse command line
    return {options.parse(argc, argv), options};
}

struct ParsedValues{
    int N_particles;
    double boxlength;
    std:: string forcefield;
    std:: string init_mode;
    double beta;
    double gamma;
    int N_iter;
    double stepsize;
    int N_meas;
    std:: string integrator;
    int threads;
    int seed;
    std:: string output_name;

};

ParsedValues processParsedValues(const cxxopts::ParseResult& result) {
    // Access parsed values
    ParsedValues values;

    values.Nparticles   =   result["N_particles"].as<int>();
    values.boxlength    =   result["boxlength"].as<double>();
    values.forcefield   =   result["forcefield"].as<std::string>();
    values.init_mode    =   result["init_mode"].as<std::string>();
    values.beta         =   result["beta"].as<double>();
    values.gamma        =   result["gamma"].as<double>();
    values.Niter        =   result["N_iter"].as<double>();
    values.stepsize     =   result["stepsize"].as<double>();
    values.N_meas       =   result["N_meas"].as<int>();
    value.integrator    =   result["integrator"].as<std::string>();
    value.threads       =   result["threads"].as<int>();
    value.seed          =   result["seed"].as<int>();
    value.output_name   =   result["output_name"].as<std::string>();

    return values;
}


#endif // ARGPARSER_H