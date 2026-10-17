#pragma once

#include "FFT.h"
#include <vector>
#include <complex>

namespace fourier{
// Struct for the result object from recover_density
    struct DensityResult {
        std::vector<double> x;        // grid points
        std::vector<double> density;  // f(x) at those points
    };

//  This function takes characteristic function (which we get from the model specific characteristic function) and the du, x_min from build_grid function return object;
//  It provides the density function and x-grid, which can then be used in the pricing function;
    DensityResult RecoverDensity(fourier::CF cf, double du, double x_min, int N = 64);
}