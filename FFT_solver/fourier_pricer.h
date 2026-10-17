#pragma once

#include "FFT.h"
#include "density_recovery.h"

namespace fourier{
// This calculates the call price by taking characteristic function, x-grid, and K,r and T as input;
    double call_price_from_density(const DensityResult& result, double K, double r, double T);
    double put_price_from_density(const DensityResult& result, double K, double r, double T);
}