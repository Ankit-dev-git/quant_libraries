#pragma once

namespace model_params{  
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

    struct MertonParams{
        double S_0,
            T,
            r,
            sigma, 
            lambda,
            mu_j,
            sigma_j;
    };

    struct BatesParams{
        HestonParams heston;
        double lambda,
            mu_j,
            sigma_j;
    };
}