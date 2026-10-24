#pragma once
#include "FFT.h"
#include <string>
#include <utility>

namespace fourier{

// This is the BS parameter struct which is passed into BS characteristic funciton builder;    
    struct BlackScholes_params{
        double S_0;
        double r;
        double sigma;
        double T;
        
    };

// This is the Heston parameter struct which is passed into Heston characteristic funciton builder;
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

    struct MertonParams{
        double S_0,
            T,
            r,
            sigma, 
            lambda,
            mu_j,
            sigma_j;
    };
// This is the struct of Characteristic funcion, which is returned from every model characteristic function builder;
    struct CharacteristicFuncResult{
        CF charac_func;
        double mu_terminal;
        double sigma_terminal;
        double c4;
        std::string model;
    };

    struct Cumulants{
        double c1,
               c2,
               c4;
    };

//  cumulant.c1 is the mean of the terminal X and cumulant.c2 is the variance of terminal X; We pass the characteristic function of X_T based on time t; X = ln(S)
//  the truncation error is about h² ≈ 1e-8, and the rounding error on c₂ is about 1e-16/h² ≈ 1e-8. Both are negligible for a window with L = 8–10 sd;
//  |h| in [1e-5, 1e-2]; negative h is valid (conjugate symmetry);
    Cumulants charac_func_cumulants(const CF& charac_func, double h = 1e-4, double h2 = 1e-2);

    CharacteristicFuncResult charac_func_BS(BlackScholes_params params);
    CharacteristicFuncResult charac_func_Heston(HestonParams params); 
    CharacteristicFuncResult charac_func_Merton(MertonParams params);

// TODO  
// 1.  
// 2. CharacteristicFuncResult charac_func_Bates(Bates params);

}