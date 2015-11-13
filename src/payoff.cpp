#ifndef PAYOFF_CPP
#define PAYOFF_CPP

#include "payoff.h"

/* ------------------- PayOff ----------------------- */

PayOff::PayOff() {}

/* ------------------- PayOffCall  ----------------------- */

PayOffCall::PayOffCall(const double &_K) : K(_K) {}

double PayOffCall::operator() (const double &S) const
{
    return std::max(S - K, 0.0);
}


/* ------------------- Class name ----------------------- */


#endif // PAYOFF_CPP
