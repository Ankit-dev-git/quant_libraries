#pragma once
#include "FFT.h"

namespace fourier{
    struct BlackScholes_params{
        double S_0;
        double r;
        double sigma;
        double T;
        
    };

    struct BSInvFourResult{
        CF charac_func;
        double mu;
    };

    BSInvFourResult charac_func_BS(BlackScholes_params params);
}