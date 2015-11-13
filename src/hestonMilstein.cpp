#ifndef HESTONMILSTEIN_CPP
#define HESTONMILSTEIN_CPP

#include "hestonMilstein.h"


/* ------------------- HestonMilstein ----------------------- */

HestonMilstein::HestonMilstein(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

HestonMilstein::~HestonMilstein() {}

void HestonMilstein::calculateVariancePath(const mat &varianceDraws, vec &variancePath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_max = std::max(variancePath(i-1), 0.0);

        variancePath(i) = variancePath(i-1) + kappa * dt * (theta - v_max) + eta * std::sqrt(v_max * dt) * varianceDraws(i-1) + (dt * (eta * eta) / 4) * (varianceDraws(i-1) * varianceDraws(i-1) - 1);
    }
}

void HestonMilstein::calculateStockPath(const mat &stockDraws, const vec &variancePath, vec &stockPath)
{
    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_max = std::max(variancePath(i-1), 0.0);

        stockPath(i) = stockPath(i-1) * std::exp( (myOption->r - 0.5 * v_max) * dt + std::sqrt(v_max * dt) * stockDraws(i-1) );
    }
}


#endif // HESTONMILSTEIN_CPP
