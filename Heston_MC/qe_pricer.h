#pragma once

#include "modelparam_structs.h"
#include <cstdint>

namespace heston_mc{
    struct HestonMCPricerResult{
        double price,
            std_error;
    };  

    heston_mc::HestonMCPricerResult heston_mc_pricer(const model_params::HestonParams& params, double K, int paths = 5000, int steps = 100, std::uint64_t seed = 1);
}
