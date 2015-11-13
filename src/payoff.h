#ifndef PAYOFF_H
#define PAYOFF_H

//#include <armadillo>
#include <algorithm>

// base payoff class
class PayOff
{
public:
    PayOff();
    virtual ~PayOff() {};

    // turn PayOff into an abstract functor
    virtual double operator () (const double &S) const = 0;
};


// derived call option class
class PayOffCall: public PayOff
{
private:
    double K;

public:
    PayOffCall(const double &_K);
    virtual ~PayOffCall() {};

    virtual double operator () (const double &S) const;
};


#endif // PAYOFF_H


