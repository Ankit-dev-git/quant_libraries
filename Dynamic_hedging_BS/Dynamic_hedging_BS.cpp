#include <iostream>
#include <thread>
#include <chrono>
#include "generate_paths.h"
#include "BS_price.h"
#include "delta_hedge.h"
#include "DataContainerCreation.h"

int main()
{
    int numPaths{ 10000 }, numSteps{ 1000 };
    double T{ 1.0 }, r{ 0.1 }, sigma{ 0.2 }, s0{ 1.0 }, K{ 0.95 };
    double dt{ T / numSteps };

    OptionType CP{OptionType::call};

    Grid stock_prices (numPaths, numSteps+1);
    stock_prices = paths::stock_paths(numPaths, numSteps, T, r, sigma, s0);
    
    Grid Pnl(numPaths, numSteps+1);

    unsigned int numThreads = std::thread::hardware_concurrency();
    std::cout << "Number of threads: " << numThreads << std::endl;
    if (numThreads == 0) numThreads = 4;

    // Creating lambda functions
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

    worker(0, numPaths);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "Single-threaded delta hedging: " << elapsed.count() << " ms\n";

    start = std::chrono::high_resolution_clock::now();

    std::vector<std::jthread> threads;
    int chunkSize = numPaths / numThreads;
    int remainder = numPaths % numThreads;
    int startPath = 0;

    for (unsigned int t = 0; t < numThreads; ++t) {
        int endPath = startPath + chunkSize + (t < static_cast<unsigned int>(remainder) ? 1 : 0);
        threads.emplace_back(worker, startPath, endPath);
        startPath = endPath;
    }

    threads.clear();

    auto computationEnd = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> computationTime = computationEnd - start;\
    std::cout << "Delta hedging took " << computationTime.count() << " ms\n";

    savedatacontainer(Pnl, "Pnl.csv");
    auto saveEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> saveTime = saveEnd - computationEnd;
    std::cout << "Total took "<< std::chrono::duration<double, std::milli>(saveEnd - start).count()<< " ms\n";

    return 0;
}
