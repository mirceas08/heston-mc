#ifndef HESTONMILSTEIN_H
#define HESTONMILSTEIN_H

#include <armadillo>
#include <cmath>
#include "hestonDiscr.h"
#include "option.h"
using namespace arma;

class HestonMilstein: public HestonDiscretization
{
public:
    HestonMilstein(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~HestonMilstein();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};


#endif // HESTONMILSTEIN_H
