#ifndef TRANSFORMEDVOLATILITY_CPP
#define TRANSFORMEDVOLATILITY_CPP

#include "transformedVolatility.h"

/* ------------------- Transformed Volatility scheme ----------------------- */

TransformedVolatility::TransformedVolatility(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

TransformedVolatility::~TransformedVolatility() {}

void TransformedVolatility::calculateVariancePath(const mat &varianceDraws, vec &volatilityPath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);
    double parenthesisTerm = (theta - (eta * eta) / (4*kappa));
    double u, gammaStar, thetaStar;

    for (int i = 1; i < vec_size; i++) {
        u = volatilityPath(i-1) + 0.5 * kappa * (parenthesisTerm / volatilityPath(i-1) - volatilityPath(i-1)) * dt;
        gammaStar = 0.5 * (volatilityPath(i-1) + u);
        thetaStar = parenthesisTerm / gammaStar;

        volatilityPath(i) = volatilityPath(i-1) + 0.5 * kappa * (thetaStar - volatilityPath(i-1)) * dt +
                            0.5 * eta * std::sqrt(dt) * varianceDraws(i-1);
    }
}

void TransformedVolatility::calculateStockPath(const mat &stockDraws, const vec &volatilityPath, vec &stockPath)
{
    double exponential;

    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        exponential = (myOption->r - 0.5 * volatilityPath(i-1)*volatilityPath(i-1)) * dt + volatilityPath(i-1) * std::sqrt(dt) * stockDraws(i-1);
        stockPath(i) = stockPath(i-1) * std::exp(exponential);
    }
}


#endif // TRANSFORMEDVOLATILITY_CPP
