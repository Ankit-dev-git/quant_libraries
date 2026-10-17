#include "characteristic_func.h"
#include <cmath>
#include <stdexcept>

// This is Black Scholes Characteristic funciton builder. It returns the characteristics function along with mean (c1), sigma (sqrt(c2)) and c4; 
// mean and sigma are calculated directly as known value without numerical approximation;
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
    result.c4 = fourier::charac_func_cumulants(result.charac_func).c4;
    return result;
}


// This is Heston Characteristic funciton builder. 
// It returns the Heston characteristics function along with c1, sqrt(c2) and c4 calculated using cumulant funciton which uses numerical approximation;
fourier::CharacteristicFuncResult fourier::charac_func_Heston (fourier::HestonParams params){
    Complex i(0.0, 1.0);
    if (params.S_0<=0 || params.v0<0 || params.T<=0 || params.gamma<=0 || std::abs(params.rho)>=1 || params.v_bar <=0 || params.kappa<=0)
        throw std::invalid_argument("charac_func_Heston: S0, gamma, kappa, v_bar and T must be greater than zero; rho must be (-1,1); v0 must not be negative!");
    if (!std::isfinite(params.S_0) || !std::isfinite(params.v0) || !std::isfinite(params.T) || !std::isfinite(params.gamma) || !std::isfinite(std::abs(params.rho)) || 
        !std::isfinite(params.v_bar) || !std::isfinite(params.kappa) || !std::isfinite(params.r))
        throw std::invalid_argument("charac_func_Heston: S0, gamma, kappa, T, v0, v_bar, r, rho must be finite!");

    fourier::CharacteristicFuncResult result;
    double logS0 = std::log(params.S_0);

    result.charac_func = [params, i, logS0] (double u) {
        Complex d = std::sqrt((params.kappa - params.rho*i*u*params.gamma)*(params.kappa - params.rho*i*u*params.gamma) + params.gamma*params.gamma*(u*u + i*u));
        Complex b1 = (params.kappa - params.rho*i*u*params.gamma + d)/(params.gamma*params.gamma);
        Complex b2 = (params.kappa - params.rho*i*u*params.gamma - d)/(params.gamma*params.gamma);
        Complex g_tilda = b2/b1;
        Complex e_minus_d_tau = std::exp(-d*params.T);
        
        Complex A_tau = (i*u*params.r)*params.T + b2*params.v_bar*params.kappa*params.T + 2*params.v_bar*params.kappa/(params.gamma*params.gamma) * 
                        std::log((1.0 - g_tilda)/(1.0 - g_tilda*e_minus_d_tau));
        Complex B_tau = b2*(1.0 - e_minus_d_tau)/(1.0 - g_tilda*e_minus_d_tau);
        Complex C = i*u;

        return std::exp(A_tau + B_tau*params.v0 + C*logS0);
    };
    result.model = "Heston";
    Cumulants cumulants = fourier::charac_func_cumulants(result.charac_func);
    result.mu_terminal  = cumulants.c1;
    result.sigma_terminal  = std::sqrt(cumulants.c2);
    result.c4 = cumulants.c4;

    return result;
}


// This funciton takes a characteristic function and calculate c1, c2 and c4 for it; It returns them in Cumulant struct object;
fourier::Cumulants fourier::charac_func_cumulants(const CF& charac_func, double h1, double h2){
    if (!std::isfinite(h1) || std::abs(h1)<1e-5 || std::abs(h1)>1e-2 || !std::isfinite(h2) || std::abs(h2)<1e-5 || std::abs(h2)>1e-2)
        throw std::invalid_argument("h should be between 1e-2:1e-5 and finite!");
    fourier::Cumulants cumulants;
    Complex log_phi_h = std::log(charac_func(h1));
    cumulants.c1 = log_phi_h.imag()/h1;
    cumulants.c2 = -2.0*log_phi_h.real()/(h1*h1);
    cumulants.c4 = 2*(std::log(charac_func(2*h2)) - 4.0*std::log(charac_func(h2))).real()/(h2*h2*h2*h2);
    return cumulants;
}
