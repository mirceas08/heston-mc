#ifndef TRANSFORMEDVOLATILITY_H
#define TRANSFORMEDVOLATILITY_H

#include <armadillo>
#include <cmath>
#include "hestonDiscr.h"
#include "option.h"
using namespace arma;

class TransformedVolatility: public HestonDiscretization
{
public:
    TransformedVolatility(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~TransformedVolatility();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};


#endif // TRANSFORMEDVOLATILITY_H
