#include <gtest/gtest.h>
#include "characteristic_func.h"
#include "FFT.h"
#include "grid_builder.h"
#include "density_recovery.h"
#include "fourier_pricer.h"
#include <stdexcept>
#include <limits>
#include <vector>
#include "BS_price.h"

TEST(MertonCf, ReducedToBS){
    fourier::MertonParams merton_params_lambda_zero {.S_0 = 100, .T = 1.0, .r = 0.00, .sigma = 0.2, .lambda = 0.0, .mu_j = -0.1, .sigma_j = 0.2};
    fourier::MertonParams merton_params_muj_sigmaj_zero {.S_0 = 100, .T = 1.0, .r = 0.00, .sigma = 0.2, .lambda = 1.0, .mu_j = 0., .sigma_j = 0.0};
    fourier::CharacteristicFuncResult result_merton = fourier::charac_func_Merton(merton_params_lambda_zero);
    fourier::CF phi_merton = result_merton.charac_func;
    fourier::CharacteristicFuncResult result_merton_1 = fourier::charac_func_Merton(merton_params_muj_sigmaj_zero);
    fourier::CF phi_merton_1 = result_merton_1.charac_func;
    fourier::BlackScholes_params bs_params {.S_0 = 100, .r = 0.0, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result_bs = fourier::charac_func_BS(bs_params);
    fourier::CF phi_bs = result_bs.charac_func;
    std::vector<double> u_vec {0.5, 1, 2, 5};
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_merton).c1, fourier::charac_func_cumulants(phi_bs).c1, 1e-10);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_merton).c2, fourier::charac_func_cumulants(phi_bs).c2, 1e-7);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_merton).c4, fourier::charac_func_cumulants(phi_bs).c4, 1e-6);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_merton_1).c1, fourier::charac_func_cumulants(phi_bs).c1, 1e-10);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_merton_1).c2, fourier::charac_func_cumulants(phi_bs).c2, 1e-7);
    EXPECT_NEAR(fourier::charac_func_cumulants(phi_merton_1).c4, fourier::charac_func_cumulants(phi_bs).c4, 1e-6);

    for (const auto& u: u_vec){
        EXPECT_NEAR(phi_merton(u).imag(), phi_bs(u).imag(), 1e-12)<< "u = "<< u;
        EXPECT_NEAR(phi_merton(u).real(), phi_bs(u).real(), 1e-12)<< "u = "<< u;
        EXPECT_NEAR(phi_merton_1(u).imag(), phi_bs(u).imag(), 1e-12)<< "u = "<< u;
        EXPECT_NEAR(phi_merton_1(u).real(), phi_bs(u).real(), 1e-12)<< "u = "<< u;
    }
}

TEST(MertonCf, AnalyticCumulants){
    fourier::MertonParams params {.S_0 = 100, .T = 1.0, .r = 0.00, .sigma = 0.2, .lambda = 1.0, .mu_j = -0.1, .sigma_j = 0.2};
    fourier::CharacteristicFuncResult result_merton = fourier::charac_func_Merton(params);
    fourier::CF phi_merton = result_merton.charac_func;
    fourier::Cumulants cumulants = fourier::charac_func_cumulants(phi_merton);
    double kappa = std::exp(params.mu_j + params.sigma_j*params.sigma_j/2) - 1;
    double c1 = std::log(params.S_0) + (params.r - params.sigma*params.sigma/2 - params.lambda*kappa)*params.T + params.lambda*params.T*params.mu_j;
    double c2 = params.sigma*params.sigma*params.T + params.lambda*params.T*(params.mu_j*params.mu_j + params.sigma_j*params.sigma_j);
    double c4 = params.lambda*params.T*(std::pow(params.mu_j,4) + 6*params.mu_j*params.mu_j*params.sigma_j*params.sigma_j + 3*std::pow(params.sigma_j,4));

    EXPECT_NEAR(cumulants.c1, c1, 1e-10);
    EXPECT_NEAR(cumulants.c2, c2, 1e-7);
    EXPECT_NEAR(cumulants.c4, c4, 1e-3*c4);
}

TEST(MertonCf, ReferencePriceSeries){
    fourier::MertonParams params {.S_0 = 100, .T = 1.0, .r = 0.05, .sigma = 0.2, .lambda = 1.0, .mu_j = -0.1, .sigma_j = 0.2};
    double K = 100.0;
    // Merton (1976) series: sum_n Poisson(lambda' T) weight * BS(S0, K, r_n, sigma_n, T)
    // Calculating Price directly using sum of poisson series and BS model and then comparing against the merton CF price
    double kappa = std::exp(params.mu_j + params.sigma_j*params.sigma_j/2) - 1;
    double lambda_p = params.lambda*(1 + kappa);
    double weight = std::exp(-lambda_p*params.T);          // n = 0 term
    double series = 0.0;
    // At a given n, merton series act like BS with adjusted mu and sigma
    for (int n = 0; n <= 50; ++n){
        if (n > 0) weight *= lambda_p*params.T/n;           // running Poisson weight, no n!
        double sigma_n = std::sqrt(params.sigma*params.sigma + n*params.sigma_j*params.sigma_j/params.T);
        double r_n = params.r - params.lambda*kappa + n*std::log(1 + kappa)/params.T;
        series += weight * bs::option_price(OptionType::call, params.S_0, K, sigma_n, 0.0, params.T, r_n);
    }

    // Fourier price: same pipeline as HestonCf.ReferencePrice
    int N = 1024;
    fourier::CharacteristicFuncResult result = fourier::charac_func_Merton(params);
    fourier::GridParams grid = fourier::build_grid(result, K, fourier::default_L, N);
    fourier::DensityResult density = fourier::RecoverDensity(result.charac_func, grid.du, grid.x_min, N);
    double price = fourier::call_price_from_density(density, K, params.r, params.T);
    // measured 4.8e-4 at N=1024 (O(dx^2)).
    EXPECT_NEAR(price, series, 1e-3);
}

TEST(MertonCf, ThrowsOnBadInputs){
    const fourier::MertonParams good {.S_0 = 100, .T = 1.0, .r = 0.05, .sigma = 0.2, .lambda = 1.0, .mu_j = -0.1, .sigma_j = 0.2};
    EXPECT_NO_THROW(fourier::charac_func_Merton(good));

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const std::vector<double> positive_bad    {0.0, -5.0, nan, inf, -inf};   // S0, T
    const std::vector<double> nonnegative_bad {-5.0, nan, inf, -inf};        // sigma, sigma_j, lambda (0 allowed)
    const std::vector<double> finite_bad      {nan, inf, -inf};              // r, mu_j (any sign allowed)

    auto expect_throws = [&](const std::vector<double>& bad_values, double fourier::MertonParams::* field, const char* name){
        for (double v : bad_values){
            fourier::MertonParams p = good;
            p.*field = v;
            EXPECT_THROW(fourier::charac_func_Merton(p), std::invalid_argument) << name << " = " << v;
        }
    };

    expect_throws(positive_bad,    &fourier::MertonParams::S_0,     "S_0");
    expect_throws(positive_bad,    &fourier::MertonParams::T,       "T");
    expect_throws(nonnegative_bad, &fourier::MertonParams::sigma,   "sigma");
    expect_throws(nonnegative_bad, &fourier::MertonParams::sigma_j, "sigma_j");
    expect_throws(nonnegative_bad, &fourier::MertonParams::lambda,  "lambda");
    expect_throws(finite_bad,      &fourier::MertonParams::r,       "r");
    expect_throws(finite_bad,      &fourier::MertonParams::mu_j,    "mu_j");
}