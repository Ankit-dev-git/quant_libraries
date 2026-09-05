#include "GridCreation.h"

Grid::Grid(int r, int c) : rows(r), cols(c), data(r* c) {}

double& Grid::operator()(int i, int j) { return data[i * cols + j]; }

double Grid::operator()(int i, int j) const { return data[i * cols + j]; }
