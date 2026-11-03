#pragma once

#include "modelparam_structs.h"

namespace heston_mc{
    inline constexpr double psi_c  = 1.5;
    struct HestonMCVarianceMoments{
    double v_next = 0,
        a = 0,
        b_sqr = 0,
        p = 0,
        beta = 0,
        delta_t = 0, 
        vt = 0;
    bool is_exp {false}; // True if exponential branch is choosen
};

// dt must be greter than zero as at dt = 0, ψ = 0, so 2/ψ = ∞ and a(√b² + Z)² = 0·∞ = NaN
    heston_mc::HestonMCVarianceMoments heston_mc_variance(double vt, double delta_t, const model_params::HestonParams& params, double u);
    double heston_mc_x_tdelta(double xt, const heston_mc::HestonMCVarianceMoments& var_moments, const model_params::HestonParams& params, double Z);
}


