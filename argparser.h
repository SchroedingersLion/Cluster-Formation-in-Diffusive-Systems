#ifndef ARGPARSER_H
#define ARGPARSER_H

#include <string>
#include <iostream>
#include <cxxopts.hpp>

// DEFAULT PARAMETERS
constexpr std:: string  _Nparticles_default {"1000"};
constexpr std:: string  _boxlength_default {"10"};
constexpr std:: string  _forcefield_default {"gauss"};
constexpr std:: string  _init_mode_default {"uniform"};
constexpr std:: string  _beta_default {"10"};
constexpr std:: string  _gamma_default {"0.1"};
constexpr std:: string  _Niter_default{"10000"};
constexpr std:: string  _stepsize_default {"0.01"};
constexpr std:: string  _Nmeas_default {"10"};
constexpr std:: string  _integrator_default {"BAOAB"};
constexpr std:: string  _threads_default {"4"};
constexpr std:: string  _seed_default {"1"};
constexpr std:: string  _output_name_default {"results.csv"};



// std:: pair <cxxopts:: ParseResult, cxxopts:: Options> parseCommandLine(int argc, char* argv[]) {
//     cxxopts::Options options("MyProgram", "Description");

//     // Define command line options
//     options.add_options()
//         ("Nparticles",  "Number of particles.",                                         cxxopts:: value <int>()->default_value("1000"))
//         ("boxlength",   "Length of edge of square simulation box.",                     cxxopts:: value <double>()->default_value("10"))
//         ("forcefield",  "Forcefield between two particles. Either 'gauss' or 'morse'.", cxxopts:: value <std:: string>()->default_value("gauss"))
//         ("init_mode",   "Initial positions of the system. Either 'uniform' or 'grid'.", cxxopts:: value <std:: string>()->default_value("uniform"))
//         ("beta",        "Inverse temperature parameter in Langevin dynamics.",          cxxopts:: value <double>()->default_value("10"))
//         ("gamma",       "Friction parameter in Langevin dynamics.",                     cxxopts:: value <double>()->default_value("0.1"))
//         ("Niter",       "Number of simulation steps.",                                  cxxopts:: value <int>()->default_value("10000"))
//         ("stepsize",    "Simulation stepsize.",                                         cxxopts:: value <double>()->default_value("0.01"))
//         ("Nmeas",       "Take a measurement any 'Nmeas' iterations.",                   cxxopts:: value <int>()->default_value("10"))
//         ("integrator",  "Integrator to be used. Either 'BAOAB' or 'UBU'.",              cxxopts:: value <std:: string>()->default_value("BAOAB"))
//         ("threads",     "Number of threads used in the force evaluation.",              cxxopts:: value <int>()->default_value("4"))
//         ("seed",        "Randomseed.",                                                  cxxopts:: value <int>()->default_value("1"))
//         ("output_name", "Name of the printed file holding the results.",                cxxopts:: value <std:: string>()->default_value("results_test.csv"))               
//         ("help",        "Print help");

//     // Parse command line
//     return {options.parse(argc, argv), options};
// }
std:: pair <cxxopts:: ParseResult, cxxopts:: Options> parseCommandLine(int argc, char* argv[]) {
    cxxopts::Options options("MyProgram", "Description");

    // Define command line options
    options.add_options()
        ("Nparticles",  "Number of particles.",                                         cxxopts:: value <int>()->default_value(_Nparticles_default))
        ("boxlength",   "Length of edge of square simulation box.",                     cxxopts:: value <double>()->default_value(_boxlength_default))
        ("forcefield",  "Forcefield between two particles. Either 'gauss' or 'morse'.", cxxopts:: value <std:: string>()->default_value(_forcefield_default))
        ("init_mode",   "Initial positions of the system. Either 'uniform' or 'grid'.", cxxopts:: value <std:: string>()->default_value(_init_mode_default))
        ("beta",        "Inverse temperature parameter in Langevin dynamics.",          cxxopts:: value <double>()->default_value(_beta_default))
        ("gamma",       "Friction parameter in Langevin dynamics.",                     cxxopts:: value <double>()->default_value(_gamma_default))
        ("Niter",       "Number of simulation steps.",                                  cxxopts:: value <int>()->default_value(_Niter_default))
        ("stepsize",    "Simulation stepsize.",                                         cxxopts:: value <double>()->default_value(_stepsize_default))
        ("Nmeas",       "Take a measurement any 'Nmeas' iterations.",                   cxxopts:: value <int>()->default_value(_Nmeas_default))
        ("integrator",  "Integrator to be used. Either 'BAOAB' or 'UBU'.",              cxxopts:: value <std:: string>()->default_value(_integrator_default))
        ("threads",     "Number of threads used in the force evaluation.",              cxxopts:: value <int>()->default_value(_threads_default))
        ("seed",        "Randomseed.",                                                  cxxopts:: value <int>()->default_value(_seed_default))
        ("output_name", "Name of the printed file holding the results.",                cxxopts:: value <std:: string>()->default_value(_output_name_default))               
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

    // if (result.count("arg1")) {
    //     arg1 = result["arg1"].as<int>();
    // } else {
    //     arg1 = default_value;
    // }
ParsedValues processParsedValues(const cxxopts:: ParseResult& result) {
    // Access parsed values
    ParsedValues values;

    values.Nparticles = result.count("Nparticles")  ?   result["Nparticles"].as<int>()         : static_cast<int>(_Nparticles_default);
    values.boxlength  = result.count("boxlength")   ?   result["boxlength"].as<double>()       : static_cast<double>(_boxlength_default);
    values.forcefield = result.count("forcefield")  ?   result["forcefield"].as<std::string>() : _forcefield_default;
    values.init_mode  = result.count("init_mode")   ?   result["init_mode"].as<std::string>()  : _init_mode_default;
    values.beta       = result.count("beta")        ?   result["beta"].as<double>()            : static_cast<double>(_beta_default);
    values.gamma      = result.count("gamma")       ?   result["gamma"].as<double>()           : static_cast<double>(_gamma_default);
    values.Niter      = result.count("Niter")       ?   result["Niter"].as<int>()              : static_cast<int>(_Niter_default);
    values.stepsize   = result.count("stepsize")    ?   result["stepsize"].as<double>()        : static_cast<double> (_stepsize_default);
    values.N_meas     = result.count("N_meas")      ?   result["Nmeas"].as<int>()              : static_cast<int>(_Nmeas_default);
    value.integrator  = result.count("integrator")  ?   result["integrator"].as<std::string>() : _integrator_default;
    value.threads     = result.count("threads")     ?   result["threads"].as<int>()            : static_cast<int>(_threads_default);
    value.seed        = result.count("seed")        ?   result["seed"].as<int>()               : static_cast<int>(_seed_default);
    value.output_name = result.count("output_name") ?   result["output_name"].as<std::string>(): _output_name_default;

    return values;
}

// ParsedValues processParsedValues(const cxxopts:: ParseResult& result) {
//     // Access parsed values
//     ParsedValues values;

//     values.Nparticles   =   result["N_particles"].as<int>();
//     values.boxlength    =   result["boxlength"].as<double>();
//     values.forcefield   =   result["forcefield"].as<std:: string>();
//     values.init_mode    =   result["init_mode"].as<std:: string>();
//     values.beta         =   result["beta"].as<double>();
//     values.gamma        =   result["gamma"].as<double>();
//     values.Niter        =   result["N_iter"].as<double>();
//     values.stepsize     =   result["stepsize"].as<double>();
//     values.N_meas       =   result["N_meas"].as<int>();
//     value.integrator    =   result["integrator"].as<std:: string>();
//     value.threads       =   result["threads"].as<int>();
//     value.seed          =   result["seed"].as<int>();
//     value.output_name   =   result["output_name"].as<std:: string>();

//     return values;
// }


#endif // ARGPARSER_H