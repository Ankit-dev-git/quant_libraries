#include <iostream>
#include <chrono>
#include "generate_paths.h"
#include "BS_price.h"
#include "delta_hedge.h"
#include "DataContainerCreation.h"
#include "parallel_for.h"

int main()
{
    int numPaths{ 10000 }, numSteps{ 1000 };
    double T{ 1.0 }, r{ 0.1 }, sigma{ 0.2 }, s0{ 1.0 }, K{ 0.95 };
    double dt{ T / numSteps };

    OptionType CP{OptionType::call};

    Grid stock_prices (numPaths, numSteps+1);
    stock_prices = paths::stock_paths(numPaths, numSteps, T, r, sigma, s0);
    
    Grid Pnl(numPaths, numSteps+1);

    // Creating a lambda funciton
    auto worker = [&](int startPath, int endPath) {
        std::vector<double> s(numSteps + 1), deltas(numSteps + 1), callPrice(numSteps + 1);
        for (int i = startPath; i < endPath; ++i) {
            for (int j = 0; j < numSteps + 1; ++j) {
                s[j] = stock_prices(i, j);
                deltas[j] = bs::option_delta(CP, s[j], K, sigma, j * dt, T, r);
                callPrice[j] = bs::option_price(CP, s[j], K, sigma, j * dt, T, r);
            }
            std::vector<double> pathPnl = delta_hedge_pnl(callPrice, deltas, s, K, r, dt);
            std::copy(pathPnl.begin(), pathPnl.end(), Pnl.data.begin() + i * Pnl.cols);
        }
    };

    auto start = std::chrono::high_resolution_clock::now();

    // Running the whole process without multi-threading
    worker(0, numPaths);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "Single-threaded delta hedging: " << elapsed.count() << " ms\n";

    // running whole process with multi-threading
    start = std::chrono::high_resolution_clock::now();

    parallel_for(0, numPaths, worker);

    auto computationEnd = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> computationTime = computationEnd - start;
    std::cout << "Multi-threaded Delta hedging took " << computationTime.count() << " ms\n";

    // Saving the final data containers in csv files for further analysis
    savedatacontainer(Pnl, "Pnl.csv");
    auto saveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> saveTime = saveEnd - computationEnd;
    std::cout << "saving took " << saveTime << " ms\n";
    std::cout << "Total took "<< std::chrono::duration<double, std::milli>(saveEnd - start).count()<< " ms\n";

    return 0;
}
