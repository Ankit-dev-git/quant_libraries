#pragma once
#include "characteristic_func.h"


namespace fourier{
    inline constexpr double default_L = 8.0;
    struct GridParams {
        double du;
        double x_min;
        double points_per_sd;
        int n_dx_xmin;
        };

// This function returns the Grid_parameters, mainly du and x_min. This parameters are then passed directly into RecoverDensity function
    GridParams build_grid(const fourier::CharacteristicFuncResult& characteristic_func_XT, double K, double L = default_L, int N = 64);
}