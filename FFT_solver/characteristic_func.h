#pragma once
#include "FFT.h"
#include "modelparam_structs.h"
#include <string>
#include <utility>


namespace fourier{

    using BlackScholes_params = model_params::BlackScholes_params;
    using HestonParams = model_params::HestonParams;
    using MertonParams = model_params::MertonParams;
    using BatesParams = model_params::BatesParams;
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
    CharacteristicFuncResult charac_func_Bates(BatesParams params);

}