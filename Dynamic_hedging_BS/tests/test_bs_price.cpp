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