#ifndef HESTONEULER_H
#define HESTONEULER_H

#include <armadillo>
#include <cmath>
#include "hestonDiscr.h"
#include "option.h"
using namespace arma;


// reflection
class HestonEulerReflection: public HestonDiscretization
{
public:
    HestonEulerReflection(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~HestonEulerReflection();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};

// partial truncation
class HestonEulerPT: public HestonDiscretization
{
public:
    HestonEulerPT(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~HestonEulerPT();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};

// full truncation
class HestonEulerFT: public HestonDiscretization
{
public:
    HestonEulerFT(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~HestonEulerFT();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};

#endif // HESTONEULER_H
