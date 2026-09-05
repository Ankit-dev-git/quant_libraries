#include "generate_paths.h"
#include "normal.h"
#include "GridCreation.h"
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>


Grid paths::stock_paths(int numPaths, int numSteps, double T, double r, double sigma, double stock_price) {
    Grid Z(numPaths, numSteps);
    Grid X(numPaths, numSteps + 1);
    Grid W(numPaths, numSteps + 1);
    Grid S(numPaths, numSteps + 1);
    std::vector<double> time_grid(numSteps + 1);

    // Z is a flattened vector of size numPaths*numSteps, not a 2D array.
    std::generate(Z.data.begin(), Z.data.end(), normal::randn);

    for (auto i = 0; i < numPaths; ++i)
        X(i, 0) = std::log(stock_price);

    double dt = T / numSteps;


    for (int j = 0; j < numSteps; ++j) {
        if (numPaths > 1) {
            double mean{ 0 }, stdev{ 0 };

            for (int i = 0; i < numPaths; ++i)
                mean += Z(i, j);
            mean /= numPaths;

            for (int i = 0; i < numPaths; ++i)
                stdev += (Z(i, j) - mean) * (Z(i, j) - mean);
            stdev = std::sqrt(stdev / numPaths);

            for (int i = 0; i < numPaths; i++)
                Z(i, j) = (Z(i, j) - mean) / stdev;
        }
        for (int i = 0; i < numPaths; ++i)
            W(i, j + 1) = W(i, j) + Z(i, j) * std::sqrt(dt);
        for (int i = 0; i < numPaths; ++i)
            X(i, j + 1) = X(i, j) + (r - sigma * sigma / 2) * dt + sigma * (W(i, j + 1) - W(i, j));
        time_grid[j + 1] = time_grid[j] + dt;
    }
    for (int i = 0; i < (numSteps + 1) * numPaths; ++i)
        S.data[i] = std::exp(X.data[i]);
    return S;
}