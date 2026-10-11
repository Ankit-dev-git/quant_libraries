#pragma once
#include "characteristic_func.h"

namespace fourier{
    
    struct GridParams {
        double du;
        double x_min;
        double points_per_sd;
        int n_dx_xmin;
        };

    GridParams build_grid(const fourier::CharacteristicFuncResult& characteristic_func_XT, double K, double L, int N = 64);
}