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
    //srand(time(NULL));
    arma_rng::set_seed_random();

    int numSims = 100000;
    int numIntervals = 1000;

    double S_0 = 100.0;    // Initial spot price
    double K = 100.0;      // Strike price
    double r = 0.0319;     // Risk-free rate
    double v_0 = 0.010201; // Initial volatility
    double T = 1.00;       // One year until expiry

    double kappa = 6.21;   // Mean-reversion rate
    double theta = 0.019;  // Long run average volatility
    double eta = 0.61;      // "Vol of vol"

    double rho = -0.7;     // Correlation of asset and volatility
    mat correlationMatrix = eye<mat>(2,2);
    correlationMatrix(0,1) = rho;
    correlationMatrix(1,0) = rho;
    mat choleskyMatrix = chol(correlationMatrix, "lower");

    PayOff* myPayoff = new PayOffCall(K);
    Option* myOption = new Option(K, r, T, myPayoff);
    HestonDiscretization* myHeston = new HestonEulerFT(myOption, kappa, theta, eta, rho);
    HestonDiscretization* myHestonReflection = new HestonEulerReflection(myOption, kappa, theta, eta, rho);
    HestonDiscretization* myHestonPT = new HestonEulerPT(myOption, kappa, theta, eta, rho);
    HestonDiscretization* myHestonMilstein = new HestonMilstein(myOption, kappa, theta, eta, rho);
    HestonDiscretization* myHestonKJ = new KahlJackel(myOption, kappa, theta, eta, rho);
    HestonDiscretization* myHestonTV = new TransformedVolatility(myOption, kappa, theta, eta, rho);

    vec stockPath = zeros<vec>(numIntervals);
    stockPath.fill(S_0);
    vec volPath = zeros<vec>(numIntervals);
    volPath.fill(v_0);

    double payoffSum;
    double optionPrice;

//    payoffSum = 0.0;
//    for (int i = 0; i < numSims; i++) {
//        mat correlatedPath = zeros<mat>(2, numIntervals);
//        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);
//
//        myHeston->calculateVariancePath(correlatedPath.row(1), volPath);
//        myHeston->calculateStockPath(correlatedPath.row(0), volPath, stockPath);
//
//        payoffSum += myOption->payoff->operator()(stockPath(numIntervals-1)) * std::exp(-r*T);
//    }
//
//    optionPrice = payoffSum / static_cast<double>(numSims);
//    cout << "Option price with full truncation: " << optionPrice << endl;
//
//    payoffSum = 0.0;
//    for (int i = 0; i < numSims; i++) {
//        mat correlatedPath = zeros<mat>(2, numIntervals);
//        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);
//
//        myHestonReflection->calculateVariancePath(correlatedPath.row(1), volPath);
//        myHestonReflection->calculateStockPath(correlatedPath.row(0), volPath, stockPath);
//
//        payoffSum += myOption->payoff->operator()(stockPath(numIntervals-1));
//    }
//
//    optionPrice = payoffSum / static_cast<double>(numSims) * std::exp(-r*T);
//    cout << "Option price with reflection: " << optionPrice << endl;
//
//    payoffSum = 0.0;
//    for (int i = 0; i < numSims; i++) {
//        mat correlatedPath = zeros<mat>(2, numIntervals);
//        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);
//
//        myHestonPT->calculateVariancePath(correlatedPath.row(1), volPath);
//        myHestonPT->calculateStockPath(correlatedPath.row(0), volPath, stockPath);
//
//        payoffSum += myOption->payoff->operator()(stockPath(numIntervals-1));
//    }
//
//    optionPrice = payoffSum / static_cast<double>(numSims) * std::exp(-r*T);
//    cout << "Option price with partial truncation: " << optionPrice << endl;
//
//    payoffSum = 0.0;
//    for (int i = 0; i < numSims; i++) {
//        mat correlatedPath = zeros<mat>(2, numIntervals);
//        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);
//
//        myHestonMilstein->calculateVariancePath(correlatedPath.row(1), volPath);
//        myHestonMilstein->calculateStockPath(correlatedPath.row(0), volPath, stockPath);
//
//        payoffSum += myOption->payoff->operator()(stockPath(numIntervals-1));
//    }
//
//    optionPrice = payoffSum / static_cast<double>(numSims) * std::exp(-r*T);
//    cout << "Option price with Milstein discretization: " << optionPrice << endl;

//    payoffSum = 0.0;
//    for (int i = 0; i < numSims; i++) {
//        mat correlatedPath = zeros<mat>(2, numIntervals);
//        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);
//
//        myHestonKJ->calculateVariancePath(correlatedPath.row(1), volPath);
//        myHestonKJ->calculateStockPath(correlatedPath, volPath, stockPath);
//
//        payoffSum += myOption->payoff->operator()(stockPath(numIntervals-1)) * std::exp(-r*T);
//    }
//
//    optionPrice = payoffSum / static_cast<double>(numSims);
//    cout << "Option price with Kahl-Jackel discretization: " << optionPrice << endl;

    payoffSum = 0.0;
    for (int i = 0; i < numSims; i++) {
        mat correlatedPath = zeros<mat>(2, numIntervals);
        generateNormalCorrelationPaths(choleskyMatrix, correlatedPath);

        myHestonTV->calculateVariancePath(correlatedPath.row(1), volPath);
        myHestonTV->calculateStockPath(correlatedPath.row(0), volPath, stockPath);

        payoffSum += myOption->payoff->operator()(stockPath(numIntervals-1)) * std::exp(-r*T);
    }

    optionPrice = payoffSum / static_cast<double>(numSims);
    cout << "Option price with Transformed Volatility scheme: " << optionPrice << endl;


    delete myPayoff;
    delete myOption;
    delete myHeston;
    delete myHestonReflection;
    delete myHestonPT;
    delete myHestonMilstein;
    delete myHestonKJ;
    delete myHestonTV;


    return 0;
}
