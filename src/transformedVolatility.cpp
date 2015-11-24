#ifndef TRANSFORMEDVOLATILITY_CPP
#define TRANSFORMEDVOLATILITY_CPP

#include "transformedVolatility.h"
#include "math.h"

/* ------------------- Transformed Volatility scheme with Central Discretization ----------------------- */

TV_centralDiscretization::TV_centralDiscretization(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

TV_centralDiscretization::~TV_centralDiscretization() {}

void TV_centralDiscretization::calculateVariancePath(const mat &varianceDraws, vec &volatilityPath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);
    double parenthesisTerm = theta - (eta * eta) / (4*kappa);
    double u, gammaStar, thetaStar;

    for (int i = 1; i < vec_size; i++) {
        u = volatilityPath(i-1) + 0.5 * kappa * (parenthesisTerm / volatilityPath(i-1) - volatilityPath(i-1)) * dt;
        gammaStar = 0.5 * (volatilityPath(i-1) + u);
        thetaStar = parenthesisTerm / gammaStar;

        volatilityPath(i) = volatilityPath(i-1) + 0.5 * kappa * (thetaStar - volatilityPath(i-1)) * dt +
                            0.5 * eta * std::sqrt(dt) * varianceDraws(i-1);
    }
}

void TV_centralDiscretization::calculateStockPath(const mat &stockDraws, const vec &volatilityPath, vec &stockPath)
{
    double exponential;
    int vec_size = stockDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);

    for (int i = 1; i < vec_size; i++) {
        exponential = (myOption->r - 0.5 * volatilityPath(i-1)*volatilityPath(i-1)) * dt + volatilityPath(i-1) * std::sqrt(dt) * stockDraws(i-1);
        stockPath(i) = stockPath(i-1) * std::exp(exponential);
    }
}


/* ------------------- Transformed Volatility scheme with Moment Matching ----------------------- */

TV_momentMatching::TV_momentMatching(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    HestonDiscretization(_myOption, _kappa, _theta, _eta, _rho) {}

TV_momentMatching::~TV_momentMatching() {}

double firstMoment(const double &theta, const double &volatility, const double &kappa, const double &dt)
{
    return theta + (volatility * volatility - theta) * std::exp(-kappa*dt);
}

double secondMoment(const double &eta, const double &kappa, const double &dt)
{
    return ( (eta*eta) / (4*kappa) ) * (1 - std::exp(-kappa*dt));
}

void TV_momentMatching::calculateVariancePath(const mat &varianceDraws, vec &volatilityPath)
{
    int vec_size = varianceDraws.size();
    double dt = myOption->T / static_cast<double>(vec_size);
    double thetaStar, beta;
    double exponential = std::exp(-0.5*kappa*dt);
    double m1, m2;

    for (int i = 1; i < vec_size; i++) {
        m1 = firstMoment(theta, volatilityPath(i-1), kappa, dt);
        m2 = secondMoment(eta, kappa, dt);

        beta = std::sqrt( std::max(m1 - m2, 0.0) );
        thetaStar = (beta - volatilityPath(i-1) * exponential) / (1 - exponential);

        volatilityPath(i) = volatilityPath(i-1) + 0.5 * kappa * (thetaStar - volatilityPath(i-1)) * dt +
                            0.5 * eta * std::sqrt(dt) * varianceDraws(i-1);
    }
}

void TV_momentMatching::calculateStockPath(const mat &stockDraws, const vec &volatilityPath, vec &stockPath)
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
