#ifndef TRANSFORMEDVOLATILITY_H
#define TRANSFORMEDVOLATILITY_H

#include <armadillo>
#include <cmath>
#include "hestonDiscr.h"
#include "option.h"
using namespace arma;

class TV_centralDiscretization: public HestonDiscretization
{
public:
    TV_centralDiscretization(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~TV_centralDiscretization();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};

class TV_momentMatching: public HestonDiscretization
{
public:
    TV_momentMatching(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~TV_momentMatching();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
    friend double firstMoment(const double &theta, const double &variance, const double &kappa, const double &dt);
    friend double secondMoment(const double &eta, const double &kappa, const double &dt);
};

#endif // TRANSFORMEDVOLATILITY_H
