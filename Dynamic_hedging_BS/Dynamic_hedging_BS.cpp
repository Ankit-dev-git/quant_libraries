#include <iostream>
#include "generate_paths.h"
#include "BS_price.h"
#include "delta_hedge.h"
#include "DataContainerCreation.h"

int main()
{
    int numPaths{ 2000 }, numSteps{ 100 };
    double T{ 1.0 }, r{ 0.1 }, sigma{ 0.2 }, s0{ 1.0 }, K{ 0.95 };
    double dt{ T / numSteps };

    OptionType CP{OptionType::call};

    Grid stock_prices (numPaths, numSteps+1);
    stock_prices = paths::stock_paths(numPaths, numSteps, T, r, sigma, s0);
    
    Grid Pnl(numPaths, numSteps+1);

    for (int i = 0; i < numPaths; ++i) {
        std::vector<double> s (numSteps + 1 ), deltas (numSteps+1), callPrice(numSteps+1);
        for (int j = 0; j < numSteps + 1; ++j) {
            s[j] = stock_prices(i, j);
            deltas[j] = bs::option_delta(CP, s[j], K, sigma, j*dt, T, r);
            callPrice[j] = bs::option_price(CP, s[j], K, sigma, j * dt, T, r);
        }
        std::vector<double> pathPnl = delta_hedge_pnl(callPrice, deltas, s, K, r, dt);
        std::copy(pathPnl.begin(), pathPnl.end(), Pnl.data.begin() + i*Pnl.cols);
   }
    
    savedatacontainer(Pnl, "Pnl.csv");

    return 0;
}
