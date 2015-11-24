#include <iostream>
#include <algorithm>
#include "option.h"
#include "payoff.h"
#include "hestonDiscr.h"
#include "hestonEuler.h"
#include "hestonMilstein.h"
#include "KahlJackel.h"
#include "transformedVolatility.h"
#include "helpers.h"

#include <armadillo>
using namespace arma;

int main(int argc, char **argv)
{
    arma_rng::set_seed_random();

    int numSims;                            // number of simulations
    int numIntervals;                       // number of time steps
    double S0;                              // initial spot price
    double K;                               // strike
    double r;                               // risk-free rate
    double v0;                              // initial variance
    double T;                               // maturity
    double kappa;                           // speed of reversion
    double theta;                           // long run average variance
    double eta;                             // vol of vol
    double rho;                             // correlation between the two processes
    std::string discretizationScheme;       // discretization scheme

    std::string dataFile = "data/data.dat";
    std::ifstream fIN(dataFile.c_str());
    std::string line;

    while (std::getline(fIN, line)) {
        std::stringstream stream(line);
        std::string variable;
        std::string value;

        stream >> variable >> value;

        if (variable == "numSims")
            numSims = static_cast<int>(atoi(value.c_str()));
        else if (variable == "timeSteps")
            numIntervals = static_cast<int>(atoi(value.c_str()));
        else if (variable == "S0")
            S0 = atof(value.c_str());
        else if (variable == "strike")
            K = atof(value.c_str());
        else if (variable == "r")
            r = atof(value.c_str());
        else if (variable == "v0")
            v0 = atof(value.c_str());
        else if (variable == "maturity")
            T = atof(value.c_str());
        else if (variable == "kappa")
            kappa = atof(value.c_str());
        else if (variable == "theta")
            theta = atof(value.c_str());
        else if (variable == "eta")
            eta = atof(value.c_str());
        else if (variable == "rho")
            rho = atof(value.c_str());
        else if (variable == "discretizationScheme")
            discretizationScheme = value;
    }

    transform(discretizationScheme.begin(), discretizationScheme.end(), discretizationScheme.begin(), ::toupper);

    /* ------------------------ Set correlation matrix and compute Cholesky ------------------------ */
    mat correlationMatrix = eye<mat>(2,2);
    correlationMatrix(0,1) = rho;
    correlationMatrix(1,0) = rho;
    mat choleskyMatrix = chol(correlationMatrix, "lower");

    /* ------------------------ Set option, payoff and discretization scheme objects ------------------------ */
    PayOff* myPayoff = new PayOffCall(K);
    Option* myOption = new Option(K, r, T, myPayoff);
    HestonDiscretization* scheme;

    if (discretizationScheme == "EULER-FT")
        scheme = new HestonEulerFT(myOption, kappa, theta, eta, rho);
    else if (discretizationScheme == "EULER-PT")
        scheme = new HestonEulerPT(myOption, kappa, theta, eta, rho);
    else if (discretizationScheme == "EULER-REFLECTION")
        scheme = new HestonEulerReflection(myOption, kappa, theta, eta, rho);
    else if (discretizationScheme == "MILSTEIN")
        scheme = new HestonMilstein(myOption, kappa, theta, eta, rho);
    else if (discretizationScheme == "KAHL-JACKEL")
        scheme = new KahlJackel(myOption, kappa, theta, eta, rho);
    else if (discretizationScheme == "TV_CD")
        scheme = new TV_centralDiscretization(myOption, kappa, theta, eta, rho);
    else if (discretizationScheme == "TV_MM")
        scheme = new TV_momentMatching(myOption, kappa, theta, eta, rho);
    else
        scheme = new HestonEulerFT(myOption, kappa, theta, eta, rho);

    vec stockPath = zeros<vec>(numIntervals+1);
    vec volPath = zeros<vec>(numIntervals+1);
    stockPath.fill(S0);
    volPath.fill(v0);

    vec payoff(numSims);
    double optionPrice;

    for (int i = 0; i < numSims; i++) {
        mat correlatedPath = zeros<mat>(2, numIntervals+1);
        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);

        if (discretizationScheme == "KAHL-JACKEL") {
            scheme->calculateVariancePath(correlatedPath.row(1), volPath);
            scheme->calculateStockPath(correlatedPath, volPath, stockPath);
        }
        else {
            scheme->calculateVariancePath(correlatedPath.row(1), volPath);
            scheme->calculateStockPath(correlatedPath.row(0), volPath, stockPath);
        }

        payoff(i) += myOption->payoff->operator()(stockPath(numIntervals)) * std::exp(-r*T);
    }

    optionPrice = sum(payoff) / static_cast<double>(numSims);
    cout << "Option price with " << discretizationScheme << ": " << optionPrice << endl;

    double sampleVar = var(payoff);
    double standardError = std::sqrt(sampleVar/numSims);

    cout << "Standard error: " << standardError << endl;

    delete myPayoff;
    delete myOption;
    delete scheme;

    return 0;
}
