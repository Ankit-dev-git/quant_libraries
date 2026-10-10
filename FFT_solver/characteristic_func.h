#pragma once
#include "FFT.h"
#include <string>

namespace fourier{
    struct BlackScholes_params{
        double S_0;
        double r;
        double sigma;
        double T;
        
    };

    struct CharacteristicFuncResult{
        CF charac_func;
        double mu_terminal;
        double sigma_terminal;
        std::string model;
    };

    CharacteristicFuncResult charac_func_BS(BlackScholes_params params);
// TODO  
// 1. CharacteristicFuncResult charac_func_Heston(Heston params); 
// 2. CharacteristicFuncResult charac_func_Merton(Merton params); 
// 3. CharacteristicFuncResult charac_func_Bates(Bates params);

}