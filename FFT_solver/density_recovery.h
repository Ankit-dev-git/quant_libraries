#pragma once

#include "FFT.h"
#include <vector>
#include <functional>
#include <complex>

// Struct for the result object from recover_density
struct DensityResult {
    std::vector<double> x;        // grid points
    std::vector<double> density;  // f(x) at those points
};

// Characteristic functions signature
using CF = std::function<Complex(double)>;

DensityResult RecoverDensity(CF cf, double du, double x_min, int N = 64);