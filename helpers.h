#ifndef HELPERS_H
#define HELPERS_H

#include "option.h"
#include "payoff.h"
#include "hestonDiscr.h"

#include <armadillo>
using namespace arma;

void generateNormalCorrelationPaths(const mat &cholesky, mat &correlatedNormalPaths)
{
    int numOfProcesses = correlatedNormalPaths.n_rows;
    int numOfPaths = correlatedNormalPaths.n_cols;

    mat uncorrelatedNormalPaths = randn<mat>(numOfProcesses, numOfPaths);

    correlatedNormalPaths = cholesky * uncorrelatedNormalPaths;
}

#endif // HELPERS_H
