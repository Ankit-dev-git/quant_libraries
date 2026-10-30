#include <gtest/gtest.h>
#include "characteristic_func.h"
#include "FFT.h"
#include "grid_builder.h"
#include "density_recovery.h"
#include "fourier_pricer.h"
#include "BS_price.h"
#include <stdexcept>
#include <limits>
#include <vector>

TEST(BatesCF, LambdaZeroHestonCf){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    fourier::BatesParams bates_params { .heston = heston_params, .lambda = 0.0, .mu_j = -0.1, .sigma_j = 0.2 };

    fourier::CharacteristicFuncResult result_bates = fourier::charac_func_Bates(bates_params);
    fourier::CF phi_bates = result_bates.charac_func;
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    fourier::CF phi_heston = result_heston.charac_func;
    std::vector<double> u_vec {0.5, 1, 2, 5};
    for (const auto& u: u_vec){
        EXPECT_NEAR(phi_bates(u).imag(), phi_heston(u).imag(), 1e-12) << "u = " << u;
        EXPECT_NEAR(phi_bates(u).real(), phi_heston(u).real(), 1e-12) << "u = " << u;
    }
}

TEST(BatesCF, AdditiveCumulants){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    fourier::BatesParams bates_params { .heston = heston_params, .lambda = 1.0, .mu_j = -0.1, .sigma_j = 0.2 };
    fourier::CharacteristicFuncResult result_bates = fourier::charac_func_Bates(bates_params);
    fourier::CharacteristicFuncResult result_heston = fourier::charac_func_Heston(heston_params);
    double kappa_j = std::exp(bates_params.mu_j + bates_params.sigma_j*bates_params.sigma_j/2) - 1;
    double mu_bates = result_heston.mu_terminal - bates_params.lambda*kappa_j*heston_params.T + bates_params.lambda*heston_params.T*bates_params.mu_j;
    double var_bates = result_heston.sigma_terminal*result_heston.sigma_terminal + bates_params.lambda*heston_params.T*(bates_params.mu_j*bates_params.mu_j + bates_params.sigma_j*bates_params.sigma_j);
    double c4_bates = result_heston.c4 + bates_params.lambda*heston_params.T*(std::pow(bates_params.mu_j, 4.0) + 6*bates_params.mu_j*bates_params.mu_j*bates_params.sigma_j*bates_params.sigma_j +
                                                                            3*std::pow(bates_params.sigma_j, 4));
    EXPECT_NEAR(result_bates.mu_terminal, mu_bates, 1e-9);
    EXPECT_NEAR(result_bates.sigma_terminal*result_bates.sigma_terminal, var_bates, 1e-6);
    EXPECT_NEAR(result_bates.c4, c4_bates, 1e-3*c4_bates);
}

TEST(BatesCF, MertonLimit){
    fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.05, .kappa = 1.5768, .gamma = 0.001, .v_bar = 0.04, .v0 = 0.04, .rho = -0.5711};
    fourier::BatesParams params { .heston = heston_params, .lambda = 1.0, .mu_j = -0.1, .sigma_j = 0.2 };
    double kappa = std::exp(params.mu_j + params.sigma_j*params.sigma_j/2) - 1;
    double lambda_p = params.lambda*(1 + kappa);
    double weight = std::exp(-lambda_p*heston_params.T);          // n = 0 term
    double series = 0.0, K = 100;
    // At a given n, merton series act like BS with adjusted mu and sigma
    for (int n = 0; n <= 50; ++n){
        if (n > 0) weight *= lambda_p*heston_params.T/n;           // running Poisson weight, no n!
        double sigma_n = std::sqrt(heston_params.v0 + n*params.sigma_j*params.sigma_j/heston_params.T);
        double r_n = heston_params.r - params.lambda*kappa + n*std::log(1 + kappa)/heston_params.T;
        series += weight * bs::option_price(OptionType::call, heston_params.S_0, K, sigma_n, 0.0, heston_params.T, r_n);
    }

    // Fourier price: same pipeline as HestonCf.ReferencePrice
    int N = 1024;
    fourier::CharacteristicFuncResult result = fourier::charac_func_Bates(params);
    fourier::GridParams grid = fourier::build_grid(result, K, fourier::default_L, N);
    fourier::DensityResult density = fourier::RecoverDensity(result.charac_func, grid.du, grid.x_min, N);
    double price = fourier::call_price_from_density(density, K, heston_params.r, heston_params.T);
    // gamma = 1e-3: limit error O(gamma) vs 1/gamma^2 cancellation; gap 4.6e-4 = Merton's FFT error at N=1024. gamma=1e-4 -> 2.1e-3, 1e-5 -> negative price..
    EXPECT_NEAR(price, series, 1e-3);
}

TEST(BatesCF, ThrowsOnBadJumpInputs){
    const fourier::HestonParams heston_params {.S_0 = 100, .T = 1, .r = 0.0, .kappa = 1.5768, .gamma = 0.5751, .v_bar = 0.0398, .v0 = 0.0175, .rho = -0.5711};
    const fourier::BatesParams good { .heston = heston_params, .lambda = 1.0, .mu_j = -0.1, .sigma_j = 0.2 };
    EXPECT_NO_THROW(fourier::charac_func_Bates(good));

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const std::vector<double> nonnegative_bad {-5.0, nan, inf, -inf};   // lambda, sigma_j (0 allowed)
    const std::vector<double> finite_bad      {nan, inf, -inf};         // mu_j (any sign allowed)

    auto expect_throws = [&](const std::vector<double>& bad_values, double fourier::BatesParams::* field, const char* name){
        for (double v : bad_values){
            fourier::BatesParams p = good;
            p.*field = v;
            EXPECT_THROW(fourier::charac_func_Bates(p), std::invalid_argument) << name << " = " << v;
        }
    };

    expect_throws(nonnegative_bad, &fourier::BatesParams::lambda,  "lambda");
    expect_throws(nonnegative_bad, &fourier::BatesParams::sigma_j, "sigma_j");
    expect_throws(finite_bad,      &fourier::BatesParams::mu_j,    "mu_j");
    // Heston fields are validated by charac_func_Heston (covered in HestonCf.ThrowsOnBadInputs)
}