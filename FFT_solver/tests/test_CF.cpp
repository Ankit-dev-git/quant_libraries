#include <gtest/gtest.h>
#include "characteristic_func.h"
#include <cmath>
#include <string>
#include <complex>
#include <vector>
#include <stdexcept>

TEST(BsCf, MuT){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    double mu_T = std::log(param.S_0) + (param.r - param.sigma*param.sigma/2)*param.T;
    EXPECT_NEAR(result.mu_terminal, mu_T, 1e-14);
}
TEST(BsCf, SigmaT){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    double sigma_T = param.sigma*std::sqrt(param.T);
    EXPECT_NEAR(result.sigma_terminal, sigma_T, 1e-14);
}
TEST(BsCf, ModelName){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result = fourier::charac_func_BS(param);
    std::string model_name = "BlackScholes";
    EXPECT_EQ(result.model, model_name);
}

TEST(BsCf, PhiModulusT1){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result1 = fourier::charac_func_BS(param);
    double u {10};
    double abs_u_10_T1 = std::abs(result1.charac_func(u));
    EXPECT_NEAR(abs_u_10_T1, std::exp(-2), 1e-14);
}
TEST(BsCf, PhiModulusT2){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 2};
    fourier::CharacteristicFuncResult result2= fourier::charac_func_BS(param);
    double u {10};
    double abs_u_10_T2 = std::abs(result2.charac_func(u));
    EXPECT_NEAR(abs_u_10_T2, std::exp(-4), 1e-14);
}
TEST(BsCf, PhiArgU){
    fourier::BlackScholes_params param {.S_0 = 100, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::CharacteristicFuncResult result1 = fourier::charac_func_BS(param);
//  u needs to be small so arg is less than Pi.
    double u = 0.5;
    double mu_T = std::log(param.S_0) + (param.r - param.sigma*param.sigma/2)*param.T;
    double arg_phi_u = std::arg(result1.charac_func(u));
    EXPECT_NEAR(arg_phi_u, u*mu_T, 1e-14);
}


TEST(BsParam, InvalidS0){
    fourier::BlackScholes_params param {.S_0 = -20.0, .r = 0.05, .sigma = 0.2, .T = 1};
    fourier::BlackScholes_params param0 {.S_0 = 0, .r = 0.05, .sigma = 0.2, .T = 1};
    EXPECT_THROW(fourier::charac_func_BS(param), std::invalid_argument);
    EXPECT_THROW(fourier::charac_func_BS(param0), std::invalid_argument);
}
TEST(BsParam, InvalidSigma){
    fourier::BlackScholes_params param {.S_0 = 100.0, .r = 0.05, .sigma = -0.2, .T = 1};
    fourier::BlackScholes_params param0 {.S_0 = 100.0, .r = 0.05, .sigma = 0, .T = 1};
    EXPECT_THROW(fourier::charac_func_BS(param), std::invalid_argument);
    EXPECT_THROW(fourier::charac_func_BS(param0), std::invalid_argument);
}
TEST(BsParam, InvalidT){
    fourier::BlackScholes_params param {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = -0.5};
    fourier::BlackScholes_params param0 {.S_0 = 100.0, .r = 0.05, .sigma = 0.2, .T = 0};
    EXPECT_THROW(fourier::charac_func_BS(param), std::invalid_argument);
    EXPECT_THROW(fourier::charac_func_BS(param0), std::invalid_argument);
}
TEST(BsParam, ValidR){
    fourier::BlackScholes_params param {.S_0 = 100.0, .r = -0.05, .sigma = 0.2, .T = 1};
    EXPECT_NO_THROW(fourier::charac_func_BS(param));
}




    