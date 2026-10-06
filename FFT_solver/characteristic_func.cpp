#include "characteristic_func.h"
#include <cmath>
#include <stdexcept>

fourier::BSInvFourResult fourier::charac_func_BS(fourier::BlackScholes_params params){
    Complex i(0.0, 1.0);
    fourier::BSInvFourResult result;
    if (params.S_0<=0 || params.sigma<=0 || params.T<=0)
        throw std::invalid_argument("charac_func_BS: S0, sigma and T must be greater than zero");

    result.mu = std::log(params.S_0) + (params.r - params.sigma*params.sigma/2)*params.T; 
    double sigma_XT = params.sigma*std::sqrt(params.T);

    result.charac_func = [=](double u){
        return std::exp(i*u*result.mu - sigma_XT*sigma_XT*u*u/2);
    };
    return result;
}
