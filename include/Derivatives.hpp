/**
 * @file Derivatives.hpp
 * @author miltos-90
 * @brief Header file for the derivatives class
 * 
 * This class is responsible for numerically evaluating the partial derivatives
 * of the cost (objective) function with respect to its input variables.
 * The derivatives are estimated according to the "Richardson's deferred approach to the limit".
 * using Ridder's algorithm [1]. The implementation has been adapted from [2] with minor
 * error fixes.
 * Sources:
 * [1] Ridders, C.J.F. (1982). Advances in Engineering Software, Vol. 4, No. 2, pp. 75 - 76.
 * [2] Isaacson, E. (1989). Numerical Recipes in C: The Art of Scientific Computing 
 *      (William H. Press, Brian P. Flannery, Saul A. Teukolsky, and William T. Vetterling); 
 * 
 * @version 1.0
 * @date 2022-04-17
 * 
 */

#pragma once
#include <vector>
#include <functional>
#include<cmath>

class Derivatives
{
    using vec      = std::vector<double>;
    using function = std::function<double (vec &x)>; // Type of the cost (objective) function

    private:

        // Data
        const static double rangeRatio;  // Percentage of input range to use as an initial stepsize
        const static double safeErr;     // Exit when current error is <safeErr> times worse than the best found to far
        const static double CON;         // Stepsize is decreased by CON at each iteration
        const static double CON2;        // Stepsize squared
        const static unsigned ntab;      // no. orders for Neville's table
        function f;                      // Function to be differentiated
        size_t noDimensions;             // Input dimensionality
        vec stepSize;                    // Initial guess for the stepsize used for the differentiation
        
        // Members

        /**
        * @brief Provides a first approximation of the derivative of the cost function around
        *    a point using a 1st order central differencing scheme.
        * 
        * @param x Point around which the derivative will be evaluated
        * @param dx Stepsize
        * @param i Dimension for which the derivative will be computed
        * @return double Approximation of the partial derivative df(x_i)/dx_i
        */
        double centralDifference(vec &x, vec &dx, unsigned i);

        /**
         * @brief Approximates the partial derivative using Ridder's algorithm.
         * 
         * @param x Point around which derivative will be evaluated
         * @param dx Stepsize
         * @return double Approximation of the partial derivative df(x_i)/dx_i
         */
        double partialDerivative(vec &x, vec &dx);

        /**
         * @brief Returns the guess for the stepsize dx (different for each dimension).
         * 
         * @param lowerBound The lower bound of the decision variables
         * @param upperBound The upper bound of the decision variables
         * @return const vec containing the stepsize to be used for each dimension
         */
        const vec getStepSize(vec &lowerBound, vec &upperBound);

    public:

        /**
         * @brief Construct a new Derivatives object
         * 
         */
        Derivatives(){};
        Derivatives(vec &lowerBound, vec &upperBound, function &func):
                        f{func},
                        noDimensions{(unsigned)lowerBound.size()}, 
                        stepSize{getStepSize(lowerBound, upperBound)}{};

        // Members

       /**
        * @brief Implements Ridder's algorithm for the evaluation of the partial derivatives
        * 
        * @param x The point around which the partial derivative will be evaluated
        * @param i The dimension for which the derivative will be evaluated
        * @return double Approximation of the partial derivative df(x_i)/dx_i
        */
        double compute(vec &x, unsigned i);
};