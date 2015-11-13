#ifndef KAHLJACKEL_H
#define KAHLJACKEL_H

#include <armadillo>
#include <cmath>
#include "hestonDiscr.h"
#include "option.h"
using namespace arma;

class KahlJackel: public HestonDiscretization
{
public:
    KahlJackel(Option* _myOption, double _kappa, double _theta, double _eta, double _rho);
    virtual ~KahlJackel();

    virtual void calculateVariancePath(const mat &draws, vec &variancePath);
    virtual void calculateStockPath(const mat &draws, const vec &variancePath, vec &stockPath);
};



#endif // KAHLJACKEL_H
