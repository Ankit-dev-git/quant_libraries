#pragma once
#include "FFT.h"

namespace fourier{
    struct BlackScholes_params{
        double S_0;
        double r;
        double sigma;
        double T;
        
    };

    CF charac_func_BS(BlackScholes_params params);
}