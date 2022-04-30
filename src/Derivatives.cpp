#include<vector>
#include<functional>
#include<cmath>
#include"Derivatives.hpp"
#include<iostream>


/*  ============================= Private Members ============================== */

// Static variables
const double   Derivatives::rangeRatio = 0.1; // Percentage of input range to use as an initial stepsize
const double   Derivatives::safeErr    = 2.0; // Exit when current error is <safeErr> times worse than the best found to far
const double   Derivatives::CON        = 1.4; // Stepsize is decreased by CON at each iteration
const double   Derivatives::CON2       = Derivatives::CON * Derivatives::CON; // Stepsize squred
const unsigned Derivatives::ntab       = 10;  // no. orders for Neville's table

// Initial stepsize (guess) for the calculation of the partial derivatives
const std::vector<double> Derivatives::getStepSize(vec &lowerBound, vec &upperBound)
{
    noDimensions = lowerBound.size();
    vec stepsize = vec(noDimensions, 0.0);
    
    for (unsigned i = 0; i < noDimensions; i++) {
        stepSize[i] = rangeRatio * std::abs(upperBound[i] - lowerBound[i]);
    }

    return stepsize;
}

// Compute derivative according to 1st order central difference scheme
double Derivatives::centralDifference(vec &x, vec &dx, unsigned dimension)
{
    vec xph(noDimensions, 0.0), xmh(noDimensions, 0.0);
    std::transform(x.begin(), x.end(), dx.begin(), std::back_inserter(xph), std::plus<double>());
    std::transform(x.begin(), x.end(), dx.begin(), std::back_inserter(xmh), std::minus<double>());
    return (f(xph) - f(xmh)) / (2 * dx[dimension]);
}

/*  ============================= Public Members ============================= */

// Compute d/dx_i[f(x_i)] for dimension <i> using Ridders' algorithm
double Derivatives::compute(vec &x, unsigned dimension)
{
    double bestError = std::numeric_limits<double>::max();
    double currentError, derivative = 0.0, fac = 1.0;
    double neville[ntab][ntab] = {{0.0}};

    // Construct dx for this dimension
    vec dx(noDimensions, 0.0); 
    dx[dimension] = stepSize[dimension];

    // Compute succesive columns in the Neville table. Will move towards smaller stepsizes and higher orders of extrapolation
    neville[0][0] = centralDifference(x, dx, dimension);
    
    for (unsigned i = 1; i <= ntab; i++) {

        // Decrease stepsize
        dx[dimension] = dx[dimension] / CON;
        neville[0][i] = centralDifference(x, dx, dimension);

        for (unsigned j = 1; j <= i; j++){ // Compute extrapolations of various orders
            fac *= CON2;
            neville[j][i] = (neville[j - 1][i] * fac - neville[j - 1][i - 1]) / (fac - 1.0);
            currentError  = std::max(std::abs(neville[j][i] - neville[j - 1][i]), std::abs(neville[j][i] - neville[j - 1][i - 1]));

            // Compare each new extrapolation to one order lower, both at the present stepsize and the previous one
            if (currentError <= bestError) {
                bestError  = currentError;
                derivative = neville[j][i];
            }
        }
        // If higher order is worse by a significant factor, quit
        if (std::abs(neville[i][i] - neville[i - 1][i - 1]) >= safeErr * bestError) { break; }
    }
    return derivative;
}
