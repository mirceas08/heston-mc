#ifndef HESTONDISC_CPP
#define HESTONDISC_CPP

#include "hestonDiscr.h"

/* ------------------- HestonDiscretization ----------------------- */

HestonDiscretization::HestonDiscretization(Option* _myOption, double _kappa, double _theta, double _eta, double _rho):
    myOption(_myOption), kappa(_kappa), theta(_theta), eta(_eta), rho(_rho) {}

HestonDiscretization::~HestonDiscretization() {}

#endif // HESTONDISC_CPP
