#ifndef TRANSFORMEDVOLATILITY_H
#define TRANSFORMEDVOLATILITY_H

#include <armadillo>
#include <cmath>
#include "hestonDiscr.h"
#include "option.h"
using namespace arma;

class TranformedVolatility: public HestonDiscretization
{
public:
    TranformedVolatility(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~TranformedVolatility();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};


#endif // TRANSFORMEDVOLATILITY_H
