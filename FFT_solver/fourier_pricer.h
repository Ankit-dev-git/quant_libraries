#pragma once

#include "FFT.h"
#include "density_recovery.h"

namespace fourier{
    double call_price_from_density(const DensityResult& result, double K, double r, double T);
}