#include "qe_pricer.h"
#include "qe_scheme.h"
#include <random>
#include <algorithm>
#include <cmath>

heston_mc::HestonMCPricerResult heston_mc::heston_mc_pricer(const model_params::HestonParams& params, double K, int paths, int steps, std::uint64_t seed){
    std::mt19937_64 rn (seed);
// U and Z are uniform and normal distribution. Use rn to draw from them.
    std::uniform_real_distribution<double> U (std::nextafter(0.0, 1.0), 1.0);
    std::normal_distribution<double> Z (0.0, 1.0);

    double delta_t = params.T/static_cast<double> (steps);
    double sum = 0, sum_sq = 0, payoff = 0;
    double xt, vt;

// For each path
    for (int i = 0; i < paths; ++i){
        xt = std::log(params.S_0), vt = params.v0;
        heston_mc::HestonMCVarianceMoments var_moments;
    // for each step:
        for (int j = 0; j < steps; ++j){
            var_moments = heston_mc::heston_mc_variance(vt, delta_t, params, U(rn));
            xt = heston_mc::heston_mc_x_tdelta(xt, var_moments, params, Z(rn));
            vt = var_moments.v_next;
        }
        payoff = std::max(std::exp(xt) - K, 0.0);
        sum += payoff;
        sum_sq += payoff*payoff;
    }

    double mean_payoff = sum/paths;
    double var_payoff = sum_sq/paths - mean_payoff*mean_payoff;
    double std_error = std::sqrt(var_payoff/paths)*std::exp(-params.r*params.T);

    heston_mc::HestonMCPricerResult result;
    result.price = mean_payoff*std::exp(-params.r*params.T);
    result.std_error = std_error;
    return result;
}