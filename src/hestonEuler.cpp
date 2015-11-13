#ifndef HESTONEULER_CPP
#define HESTONEULER_CPP

#include "hestonEuler.h"


/* ------------------- HestonEulerReflection ----------------------- */

HestonEulerReflection::HestonEulerReflection(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

HestonEulerReflection::~HestonEulerReflection() {}

void HestonEulerReflection::calculateVariancePath(const mat &varianceDraws, vec &variancePath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_abs = std::abs(variancePath(i-1));

        variancePath(i) = v_abs + kappa * dt * (theta - v_abs) + eta * std::sqrt(v_abs * dt) * varianceDraws(i-1);
    }
}

void HestonEulerReflection::calculateStockPath(const mat &stockDraws, const vec &variancePath, vec &stockPath)
{
    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_abs = std::abs(variancePath(i-1));

        stockPath(i) = stockPath(i-1) * std::exp( (myOption->r - 0.5 * v_abs) * dt + std::sqrt(v_abs * dt) * stockDraws(i-1) );
    }
}


/* ------------------- HestonEulerPT ----------------------- */

HestonEulerPT::HestonEulerPT(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

HestonEulerPT::~HestonEulerPT() {}

void HestonEulerPT::calculateVariancePath(const mat &varianceDraws, vec &variancePath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_max = std::max(variancePath(i-1), 0.0);

        variancePath(i) = variancePath(i-1) + kappa * dt * (theta - variancePath(i-1)) + eta * std::sqrt(v_max * dt) * varianceDraws(i-1);
    }
}

void HestonEulerPT::calculateStockPath(const mat &stockDraws, const vec &variancePath, vec &stockPath)
{
    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_max = std::max(variancePath(i-1), 0.0);

        stockPath(i) = stockPath(i-1) * std::exp( (myOption->r - 0.5 * v_max) * dt + std::sqrt(v_max * dt) * stockDraws(i-1) );
    }
}


/* ------------------- HestonEulerFT ----------------------- */

HestonEulerFT::HestonEulerFT(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

HestonEulerFT::~HestonEulerFT() {}

void HestonEulerFT::calculateVariancePath(const mat &varianceDraws, vec &variancePath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_max = std::max(variancePath(i-1), 0.0);

        variancePath(i) = variancePath(i-1) + kappa * dt * (theta - v_max) + eta * std::sqrt(v_max * dt) * varianceDraws(i-1);
    }
}

void HestonEulerFT::calculateStockPath(const mat &stockDraws, const vec &variancePath, vec &stockPath)
{
    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        double v_max = std::max(variancePath(i-1), 0.0);

        stockPath(i) = stockPath(i-1) * std::exp( (myOption->r - 0.5 * v_max) * dt + std::sqrt(v_max * dt) * stockDraws(i-1) );
    }
}


#endif // HESTONEULER_CPP
