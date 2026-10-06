#include "characteristic_func.h"
#include <cmath>
#include <stdexcept>

fourier::CF fourier::charac_func_BS(fourier::BlackScholes_params params){
    Complex i(0.0, 1.0);
    if (params.S_0<=0 || params.sigma<=0 || params.T<=0)
        throw std::invalid_argument("charac_func_BS: S0, sigma and T must be greater than zero");

    double mu_XT = std::log(params.S_0) + (params.r - params.sigma*params.sigma/2)*params.T; 
    double sigma_XT = params.sigma*std::sqrt(params.T);

    CF char_func = [=](double u){
        return std::exp(i*u*mu_XT - sigma_XT*sigma_XT*u*u/2);
    };
    return char_func;
}
