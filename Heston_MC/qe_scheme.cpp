#include "qe_scheme.h"
#include "normal.h"
#include <cmath>

heston_mc::HestonMCVarianceMoments heston_mc::heston_mc_variance(double vt, double delta_t, const model_params::HestonParams& params, double u){
    heston_mc::HestonMCVarianceMoments var_moments;
    double kappa = params.kappa, v_bar = params.v_bar, gamma = params.gamma;
    double e = std::exp(-delta_t*kappa);
    double m = e*(vt - v_bar) + v_bar;
    double var = gamma*gamma/kappa*(v_bar/2 * (1 - e)*(1-e) + vt*e*(1-e));
//for a given set of vt, dt, and calibrated parameters - vbar, gamma, kappa: psi is deterministic;
    double psi = var/(m*m); 
    var_moments.delta_t = delta_t;
    var_moments.vt = vt;

    if (psi< heston_mc::psi_c){
        var_moments.b_sqr = 2/psi - 1 + std::sqrt(2/psi*(2/psi -1));
        var_moments.a = m/(1+var_moments.b_sqr);
        double Z = normal::norm_inv_cdf(u);
        var_moments.v_next = var_moments.a*(std::sqrt(var_moments.b_sqr) + Z)*(std::sqrt(var_moments.b_sqr) + Z);
        return var_moments;
    }

    else{
        var_moments.is_exp = true;
        var_moments.p = (psi - 1)/(psi + 1);
        var_moments.beta = (1-var_moments.p)/m;
        if (u<var_moments.p){
            var_moments.v_next = 0;
            return var_moments;
        }
        else{
            var_moments.v_next = std::log((1-var_moments.p)/(1-u))/var_moments.beta;
            return var_moments;
        }
    }
}

double heston_mc::heston_mc_x_tdelta(double xt, const heston_mc::HestonMCVarianceMoments& var_moments, const model_params::HestonParams& params, double Z){
    double k0, k1, k2, k3, k4;
    
    k1 = var_moments.delta_t/2 *(params.rho*params.kappa/params.gamma - 0.5) - params.rho/params.gamma;
    k2 = var_moments.delta_t/2 *(params.rho*params.kappa/params.gamma - 0.5) + params.rho/params.gamma;
    k3 = k4 = var_moments.delta_t/2*(1-params.rho*params.rho);
    double mgf, A = k2 + k4/2;

    if (var_moments.is_exp){
        if (A<var_moments.beta){
            mgf = var_moments.p + (1-var_moments.p)*var_moments.beta/(var_moments.beta - A);
            k0 = params.r*var_moments.delta_t - (k1+k3/2)*var_moments.vt - std::log(mgf);        
        }
        else
            k0 = (params.r - params.rho*params.kappa*params.v_bar/params.gamma)*var_moments.delta_t;
    }
    else{
        if (2*var_moments.a*A < 1){
            mgf = std::exp((A*var_moments.b_sqr*var_moments.a)/(1 - 2*A*var_moments.a))/std::sqrt((1 - 2*A*var_moments.a));
            k0 = params.r*var_moments.delta_t - (k1+k3/2)*var_moments.vt - std::log(mgf);        
        }
        else
            k0 = (params.r - params.rho*params.kappa*params.v_bar/params.gamma)*var_moments.delta_t;
    }   
    
    return xt + k0 + k1*var_moments.vt + k2*var_moments.v_next + std::sqrt(k3*var_moments.vt + k4*var_moments.v_next)*Z;
}
