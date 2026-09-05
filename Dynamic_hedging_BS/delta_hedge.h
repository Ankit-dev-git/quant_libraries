#pragma once
#include <vector>

std::vector<double> delta_hedge_pnl(const std::vector<double>& call_prices, const std::vector<double>& deltas, const std::vector<double>& stock_prices, double k, double r, double dt);
