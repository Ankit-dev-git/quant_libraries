#include <gtest/gtest.h>
#include <cmath>
#include "BS_price.h"

// Test 1: a call with a known textbook price
TEST(BlackScholes, CallKnownValue) {
    double price = bs::option_price(OptionType::call, 100.0, 100.0, 0.20, 0.0, 1.0, 0.05);
    EXPECT_NEAR(price, 10.4506, 1e-4);
}

// Test 2: put-call parity  C - P = S - K * exp(-r (T - t))
TEST(BlackScholes, PutCallParity) {
    double S = 100.0, K = 95.0, sigma = 0.25, t = 0.0, T = 0.5, r = 0.03;
    double call = bs::option_price(OptionType::call, S, K, sigma, t, T, r);
    double put  = bs::option_price(OptionType::put,  S, K, sigma, t, T, r);
    EXPECT_NEAR(call - put, S - K * std::exp(-r * (T - t)), 1e-10);
}

TEST(BlackScholes, AtExpiry){
    double S_ITM = 100.0, S_OTM = 90.0, S_ATM=95.0, K = 95.0, sigma = 0.25, t = 0.5, T = 0.5, r = 0.03;
    double call_ITM = bs::option_price(OptionType::call, S_ITM, K, sigma, t, T, r);
    double call_OTM = bs::option_price(OptionType::call, S_OTM, K, sigma, t, T, r);
    double call_ATM = bs::option_price(OptionType::call, S_ATM, K, sigma, t, T, r);
    double put_OTM = bs::option_price(OptionType::put, S_ITM, K, sigma, t, T, r);
    double put_ITM = bs::option_price(OptionType::put, S_OTM, K, sigma, t, T, r);
    double put_ATM = bs::option_price(OptionType::put, S_ATM, K, sigma, t, T, r);
    double call_ITM_delta = bs::option_delta(OptionType::call, S_ITM, K, sigma, t, T, r);
    double call_OTM_delta = bs::option_delta(OptionType::call, S_OTM, K, sigma, t, T, r);
    double call_ATM_delta = bs::option_delta(OptionType::call, S_ATM, K, sigma, t, T, r);
    double put_ITM_delta  = bs::option_delta(OptionType::put,  S_OTM, K, sigma, t, T, r);
    double put_OTM_delta  = bs::option_delta(OptionType::put,  S_ITM, K, sigma, t, T, r);
    double put_ATM_delta  = bs::option_delta(OptionType::put,  S_ATM, K, sigma, t, T, r);
    EXPECT_NEAR(call_ITM, 5.0, 1e-10);
    EXPECT_NEAR(call_OTM, 0.0, 1e-10);
    EXPECT_NEAR(call_ATM, 0.0, 1e-10);
    EXPECT_NEAR(put_ITM, 5.0, 1e-10);
    EXPECT_NEAR(put_ATM, 0.0, 1e-10);
    EXPECT_NEAR(put_OTM, 0.0, 1e-10);
    EXPECT_NEAR(call_ITM_delta, 1, 1e-10);
    EXPECT_NEAR(call_OTM_delta, 0, 1e-10);
    EXPECT_NEAR(call_ATM_delta, 0, 1e-10);
    EXPECT_NEAR(put_ITM_delta, -1, 1e-10);
    EXPECT_NEAR(put_OTM_delta, 0, 1e-10);
    EXPECT_NEAR(put_ATM_delta, 0, 1e-10);
}