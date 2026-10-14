#pragma once
#include "FFT.h"
#include <string>
#include <utility>

namespace fourier{
    struct BlackScholes_params{
        double S_0;
        double r;
        double sigma;
        double T;
        
    };

    struct HestonParams{
        double S_0,             // Initial price
        T,                      // Maturity time
        r,                      // interest rate
        kappa,                  // speed of reversion for variance
        gamma,                  // vol of vol
        v_bar,                  // long-run mean of variance
        v0,                     // initial variance
        rho;                    // correlation between X-brownian and V-brownian series
    };

    struct CharacteristicFuncResult{
        CF charac_func;
        double mu_terminal;
        double sigma_terminal;
        std::string model;
    };

//  pair.first is the mean of the terminal X and pair.second is the variance of terminal X; We pass the characteristic function of X_T based on time t; X = ln(S)
//  the truncation error is about h² ≈ 1e-6, and the rounding error on c₂ is about 1e-16/h² ≈ 1e-10. Both are negligible for a window with L = 8–10 sd;
//  |h| in [1e-5, 1e-2]; negative h is valid (conjugate symmetry);
    std::pair<double, double> charac_func_cumulants(const CF& charac_func, double h = 1e-3);

    CharacteristicFuncResult charac_func_BS(BlackScholes_params params);
    CharacteristicFuncResult charac_func_Heston(HestonParams params); 

// TODO  
// 1. CharacteristicFuncResult charac_func_Merton(Merton params); 
// 2. CharacteristicFuncResult charac_func_Bates(Bates params);

}