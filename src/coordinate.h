#ifndef COORDINATE_H
#define COORDINATE_H

#include <array>
#include <iostream>
#include <stdexcept>

template <size_t DIMENSION>
class coordinate {

    public:
        
        // Default constructor.
        coordinate() 
        {
                coords.fill(0.0); // Initialize all elements to 0.0.
        }

        // Access operator for reading/writing within the runtime dimension.
        double& operator[](size_t index) {
            return coords[index];
        }

        const double& operator[](size_t index) const {
            return coords[index];
        }

    private:
        std::array<double, DIMENSION> coords; // Fixed-size array.

};

#endif // COORDINATE_H