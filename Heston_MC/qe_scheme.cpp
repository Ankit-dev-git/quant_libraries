#include "qe_scheme.h"
#include "normal.h"
#include <cmath>

double heston_mc::heston_mc_variance(double vt, double dt, const model_params::HestonParams& params, double u){
    double kappa = params.kappa, v_bar = params.v_bar, gamma = params.gamma;
    double e = std::exp(-dt*kappa);
    double m = e*(vt - v_bar) + v_bar;
    double var = gamma*gamma/kappa*(v_bar/2 * (1 - e)*(1-e) + vt*e*(1-e));
//for a given set of vt, dt, and calibrated parameters - vbar, gamma, kappa: psi is deterministic;
    double psi = var/(m*m); 
    
    if (psi< heston_mc::psi_c){
        double b_sqr = 2/psi - 1 + std::sqrt(2/psi*(2/psi -1)), a = m/(1+b_sqr);
        double Z = normal::norm_inv_cdf(u);
        return a*(std::sqrt(b_sqr) + Z)*(std::sqrt(b_sqr) + Z);
    }

    else{
        double p = (psi - 1)/(psi + 1);
        if (u<p)
            return 0.0;
        else
            return std::log((1-p)/(1-u))/((1-p)/m);
    }
}