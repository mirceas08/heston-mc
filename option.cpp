#ifndef OPTION_CPP
#define OPTION_CPP

#include "option.h"

Option::Option(double _K, double _r, double _T, PayOff* _payoff):
    K(_K), r(_r), T(_T), payoff(_payoff) {}

Option::~Option() {}

#endif // OPTION_CPP
