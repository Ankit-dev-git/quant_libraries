#pragma once
#include<vector>

struct Grid {
    int rows, cols;
    std::vector<double> data;

    Grid(int r, int c);
    double& operator()(int i, int j);
    double  operator()(int i, int j) const;
};


