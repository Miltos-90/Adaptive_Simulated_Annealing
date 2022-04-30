/**
 * @file TemperatureScheduler.hpp
 * @author miltos-90
 * @brief Header file for the Temperature Scheduler. It implements the cooling schedule of the adaptive simulated annealing algorithm.
 * @version 0.1
 * @date 2022-04-17
 * 
 */

#pragma once
#include"Derivatives.hpp"
#include"utils.hpp"
#include<algorithm>
#include<vector>
#include<cmath>
#include<random>

class TemperatureScheduler{
    using vec  = std::vector<double>;
    using bvec = std::vector<bool>;
    using strc = utils::pointStruct<double>;     // Struct containing a point in the N-dimensional input space and it associated cost value.
    using fcn  = std::function<double (vec &x)>; // Type of the cost (objective) function

    private:
        Derivatives diff;                 // Object for the computation of derivatives
        size_t   noDimensions;            // Input dimensionality of the function that is being minimised
        double   invDimensions;           // = 1 / noDimensions
        double   m, n;                    // Parameters used to compute the temperatureParameterRatio
        unsigned numberCostSamples;       // Number of cost samples used to estimate initial temperature
        double   costParameterScaleRatio; // Ratio of cost parameter temperature annealing scales
        fcn      costFunction;            // Function to be minimised
        bvec     integerVars;             // Integer variables

        // Cost temperature-related variables
        double   currentCostTemperature = 0.0, initialCostTemperature = 0.0, costParameterRatio;  
        unsigned costIterationCounter = 0;
        
        // Variable temperature-related variables
        double temperatureParameterRatio;
        vec currentTemperatures, initialTemperatures, iterationCounters, sensitivities;
        
        // Members

       /**
        * @brief Computes the initial cost temperatures by calling the objective function <noRepeats> times with randomly sampled points
        * 
        * @param noRepeats number of times to call the cost function
        * @param lowerBound the lower bound of the input space
        * @param upperBound the upper bound of the input space
        * @param rEngine seeder for the random selection of points to sample the cost function
        */
        void getInitialTemperature(unsigned noRepeats, vec &lowerBound, vec &upperBound, std::minstd_rand &rEngine);

       /**
        * @brief Reanneals the cost temperature
        * 
        * @param bestCost Cost (objective function) value of the best point found so far
        * @param currentCost Cost (objective function) value of the current (last accepted) point
        */
        void reAnnealInitialCostTemperature(double bestCost, double currentCost);

       /**
        * @brief Reanneals the counter used for the cost temperature evaluation
        * 
        * @param bestCost Cost (objective function) value of the best point found so far
        * @param currentCost Cost (objective function) value of the current (last accepted) point
        */
        void   reAnnealCostIterationCounter(double bestCost, double currentCost);

       /**
        * @brief Reanneals the temperature for each input dimension
        * 
        * @param point Current (last accepted) point
        */
        void reAnnealTemperatures(vec &point);

    public:

        /**
         * @brief Construct a new Temperature Scheduler object
         * 
         */
        TemperatureScheduler(){};
        TemperatureScheduler(double ratioScale, 
                             double annealScale, 
                             double costScaleRatio,
                             unsigned noCostSamples):
                                m{-std::log10(ratioScale)},
                                n{std::log10(annealScale)},
                                numberCostSamples{noCostSamples},
                                costParameterScaleRatio{costScaleRatio}{};
        
        // Members

       /**
        * @brief Initialises required properties before use.
        * 
        * @param integerVariables Vector indicating if a specific dimension corresponds to an integer input space.
        * If yes, several re-annealing related operations are bypassed.
        * @param lowerBound The lower bound of the input space
        * @param upperBound The upper bound of the input space
        * @param func Cost (objective) function
        * @param rEngine Seeder for the random distributions
        */
        void setProperties(bvec &integerVariables, vec &lowerBound, vec &upperBound, fcn &func, std::minstd_rand &rEngine);

       /**
        * @brief Implements the re-annealing related operations. Wrapper of all related functions.
        * 
        * @param bestCost Cost (objective function) value of the best point found so far
        * @param currentPoint The current (last accepted) point.
        */
        void reAnneal(double bestCost, strc &currentPoint);

        /**
         * @brief Get the Cost Temperature parameter
         * 
         * @return double* cost temperature
         */
        double* getCostTemperature(void);

        /**
         * @brief Get the Temperatures parametr
         * 
         * @return vec* the variable-related temperatures
         */
        vec* getTemperatures(void);

        /**
         * @brief Increments all the required counters and temperatures at each iteration
         * 
         */
        void step(void);

        /**
         * @brief Indicates if the input dimension corresponds to an integer variable
         * 
         * @param dimension to be checked
         * @return true if it is an integer dimension
         * @return false otherwise
         */
        bool isInteger(unsigned dimension);

};