#include "characteristic_func.h"
#include <cmath>
#include <stdexcept>

fourier::CharacteristicFuncResult fourier::charac_func_BS(fourier::BlackScholes_params params){
    Complex i(0.0, 1.0);
    fourier::CharacteristicFuncResult result;
    if (params.S_0<=0 || params.sigma<=0 || params.T<=0)
        throw std::invalid_argument("charac_func_BS: S0, sigma and T must be greater than zero");

    result.mu_terminal = std::log(params.S_0) + (params.r - params.sigma*params.sigma/2)*params.T; 
    result.sigma_terminal = params.sigma*std::sqrt(params.T);
    result.model = "BlackScholes";
    result.charac_func = [mu_t=result.mu_terminal, sigma_t = result.sigma_terminal, i](double u){
        return std::exp(i*u*mu_t- sigma_t*sigma_t*u*u/2);
    };
    return result;
}
