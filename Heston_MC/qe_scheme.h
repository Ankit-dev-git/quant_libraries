#pragma once

#include "modelparam_structs.h"

namespace heston_mc{
// dt must be greter than zero as at dt = 0, ψ = 0, so 2/ψ = ∞ and a(√b² + Z)² = 0·∞ = NaN
    double heston_mc_variance(double vt, double dt, const model_params::HestonParams& params, double u);
    inline constexpr double psi_c  = 1.5;
}


