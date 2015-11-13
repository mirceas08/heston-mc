#ifndef OPTION_H
#define OPTION_H

#include "payoff.h"

class Option
{
public:
    PayOff* payoff;
    double K;
    double r;
    double T;

    Option(double _K, double _r, double _T, PayOff* _payoff);
    virtual ~Option();
};

#endif // OPTION_H
