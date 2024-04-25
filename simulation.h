#ifndef SIMULATION_H
#define SIMULATION_H

#define _USE_MATH_DEFINES

#include "model.h"
#include "measurement.h"

#include <omp.h>
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>



class simulation {

    public: 

        // CONSTRUCTOR.
        simulation (const double beta=300, const int threads=4, const std:: string integrator="BAOAB", const std:: string init_mode="uniform", const int randomseed=1 )
            : beta {beta}, THREADS {threads}, integrator {integrator}, init_mode {init_mode}, randomseed {randomseed}
            {

                // Specify force field.
                if (integrator=="BAOAB") integrator_step = &BAOAB_step;
                else if (integrator=="UBU") integrator_step = &UBU_step;
                else throw std:: invalid_argument( "Invalid integrator in simulation construction. Allowed are 'BAOAB' and 'UBU'." );

            }; 



    private:

        const double beta;
        const int THREADS;
        const std:: string integrator;  
        const std:: string init_mode;
        const int randomseed;
        void (*integrator_step) ()


};

#endif // SIMULATION_H