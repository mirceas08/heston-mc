#ifndef KAHLJACKEL_CPP
#define KAHLJACKEL_CPP

#include "KahlJackel.h"


/* ------------------- Kahl-Jackel implicit Milstein scheme ----------------------- */

KahlJackel::KahlJackel(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

KahlJackel::~KahlJackel() {}

void KahlJackel::calculateVariancePath(const mat &varianceDraws, vec &variancePath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);
    double numerator;
    double denominator = 1 + kappa * dt;
    double kappaThetaDt = kappa * theta * dt;

    for (int i = 1; i < vec_size; i++) {
        if (variancePath(i-1) < 0.0) {
            double v_max = std::max(variancePath(i-1), 0.0);

            variancePath(i) = variancePath(i-1) + kappa * dt * (theta - v_max) + eta * std::sqrt(v_max * dt) * varianceDraws(i-1);
        }
        else {
            numerator = variancePath(i-1) + kappaThetaDt + eta * std::sqrt(variancePath(i-1)*dt) * varianceDraws(i-1) + 0.25 * eta * eta * dt * (varianceDraws(i-1)*varianceDraws(i-1) - 1);
            variancePath(i) = numerator / denominator;
        }
    }
}

void KahlJackel::calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath)
{
    mat varianceDraws = draws.row(1);
    mat stockDraws = draws.row(0);

    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);
    double exponential;
    double halfRhoEtaDt = -0.25 * rho * eta * dt;

    for (int i = 1; i < vec_size; i++) {
        double vCurrent_max = std::max(variancePath(i), 0.0);
        double vLast_max = std::max(variancePath(i-1), 0.0);

        exponential = (myOption->r -0.25 * (vLast_max + vCurrent_max)) * dt + rho * std::sqrt(vLast_max*dt) * varianceDraws(i-1) +
                                0.5 * sqrt(dt) * (std::sqrt(vLast_max) + std::sqrt(vCurrent_max)) * (stockDraws(i-1) - rho * varianceDraws(i-1)) +
                                halfRhoEtaDt * (varianceDraws(i-1)*varianceDraws(i-1) - 1);
        stockPath(i) = stockPath(i-1) * std::exp(exponential);
    }

//    for (int i = 1; i < vec_size; i++) {
//        double v_max = std::max(variancePath(i-1), 0.0);
//
//        stockPath(i) = stockPath(i-1) * std::exp( (myOption->r - 0.5 * v_max) * dt + std::sqrt(v_max * dt) * stockDraws(i-1) );
//    }
}

#endif // KAHLJACKEL_CPP
