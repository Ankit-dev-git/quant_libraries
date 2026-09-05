#pragma once
#include"GridCreation.h"
#include <vector>


namespace paths {
    Grid stock_paths(int numPaths, int numSteps, double T, double r, double sigma, double stock_price);
}

