#include "delta_hedge.h"
#include "BS_price.h"
#include <cmath>
#include <algorithm>

std::vector<double> delta_hedge_pnl(const std::vector<double>& call_prices,const std::vector<double>& deltas,const std::vector<double>& stock_prices, double k, double r, double dt)
{
    int numSteps = stock_prices.size();
    std::vector<double> Pnl(numSteps);
    for (int i = 0; i < numSteps; ++i) {
        if (i == 0)
            Pnl[i] = call_prices[i] - deltas[i] * stock_prices[i];
        else 
            Pnl[i] = Pnl[i - 1] * (std::exp(r * dt)) - (deltas[i] - deltas[i - 1]) * stock_prices[i];
    }
    Pnl[numSteps - 1] = Pnl[numSteps-1] - std::max(stock_prices[numSteps - 1] - k, 0.0) +deltas[numSteps - 1] * stock_prices[numSteps - 1];
    return Pnl;
}
