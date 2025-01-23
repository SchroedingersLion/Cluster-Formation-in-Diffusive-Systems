#ifndef COORDINATE_H
#define COORDINATE_H

#include <array>
#include <iostream>
#include <stdexcept>

constexpr size_t MAXDIM = 3; // Maximum allowed dimension.

class coordinate {

    public:
        
        // Default constructor.
        coordinate() 
            : dimension {0} {
                coords.fill(0.0); // Initialize all elements to 0.0.
        }

        coordinate(size_t dim)
            : dimension {dim} {
                coords.fill(0.0);
        }

        // Method to set the runtime dimension.
        void set_dimension(size_t dim) {
            if (dim > MAXDIM) {
                throw std::out_of_range("Dimension exceeds MAXDIM.");
            }
            dimension = dim;
        }

        // Access operator for reading/writing within the runtime dimension.
        double& operator[](size_t index) {
            // if (index >= dimension) {
            //     throw std::out_of_range("Index exceeds runtime dimension.");
            // }
            return coords[index];
        }

        const double& operator[](size_t index) const {
            // if (index >= dimension) {
            //     throw std::out_of_range("Index exceeds runtime dimension.");
            // }
            return coords[index];
        }

    private:
        std::array<double, MAXDIM> coords; // Fixed-size array.
        size_t dimension;                  // Actual dimension.

        // Getters
        size_t size() const { return dimension; }

};

#endif // COORDINATE_H