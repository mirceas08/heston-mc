#ifndef HESTONDISC_H
#define HESTONDISC_H

#include <armadillo>
#include <cmath>
#include "option.h"
using namespace arma;

// base class
class HestonDiscretization
{
protected:
    Option* myOption;
    double kappa;
    double theta;
    double eta;
    double rho;
public:
    HestonDiscretization(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~HestonDiscretization();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath) = 0;
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath) = 0;
};


#endif // HESTONDISC_H
