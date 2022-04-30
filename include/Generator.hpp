/**
 * @file Generator.hpp
 * @author miltos-90
 * @brief Header file for the Generator class. 
 * 
 * This class generates new (candidate) points at each iteration of the 
 * ASA algorithm randomly sampled from a parameterised Cauchy distribution.
 * 
 * @version 1.0
 * @date 2022-04-17
 * 
 */

#pragma once
#include<vector>
#include<random>
#include"utils.hpp"

class Generator
{
    using vec  = std::vector<double>;
    using bvec = std::vector<bool>;
    using strc = utils::pointStruct<double>;                // Struct containing a point in the N-dimensional input space and it associated cost value.
    using fcn  = std::function<double (vec &x)>;            // Type of the cost (objective) function
    using bfcn = std::function<std::vector<bool> (vec &x)>; // Type of the contraint function

    private:
        unsigned frequencyModulus;       // Frequency for which re-annealing is performed
        unsigned maxStates;              // Max points to be generated (stopping criteria)
        vec lowerBound, upperBound;      // Lower and upper bounds for the decision variables
        size_t noDimensions;             // Dimensionality of the problem
        unsigned totalStates    = 0;     // Counter for total number of points generated in total
        unsigned reannealStates = 0;     // Counter for number of points generated since last re-anneal.
        fcn costFunction;                // Function to be minimised
        bfcn constraintFunction;         // Constraints function
        bvec integerVars;                // Vector containing the integer variables
        std::minstd_rand Seed;           // Seeder for random distribution

        // Members

       /**
        * @brief Checks if a given point lies outside the search space of the problem.
        * 
        * @param point N-dimensional point to be checked.
        * @return true if point lies outside the search space 
        * @return false otherwise
        */
        bool isOutOfBounds(vec &point);

        /* violatesConstraints(): 
            Inputs:
                point: 
            Outputs:
                Boolean value indicating if point violates the contraints (=1) or not (=0)
        */

       /**
        * @brief Checks if a given point violates the constraints of the problem.
        * 
        * @param point N-dimensional point to be checked.
        * @return true if point violates the contraints
        * @return false otherwise
        */
        bool violatesConstraints(vec &point);

        /**
         * @brief Checks if the input dimension corresponds to an integer variable
         * 
         * @param dimension the dimension to be checked
         * @return true if it corresponds to an integer variable
         * @return false otherwise
         */
        bool isInteger(unsigned dimension);

    public:
        
        /**
         * @brief Construct a new Generator object
         * 
         */
        Generator(){};
        Generator(unsigned int modulus, unsigned int states):
            frequencyModulus{modulus}, maxStates{states}{};

        // Members
        
        /**
         * @brief Sets uninitialised properties of the class prior to its use
         * 
         * @param intVariables vector indicating integer variables
         * @param lbound the lower bound of the search space
         * @param ubound the upper bound of the search space
         * @param func the cost (objective) function to be minimised
         * @param constFunction the constraints function
         * @param rEngine seeder for the random distributions
         */
        void setProperties(bvec &intVariables, vec &lbound, vec &ubound, fcn &func, bfcn &constFunction, std::minstd_rand &rEngine);
        
       /**
        * @brief Generates a new (candidate) point sampled from a parameterised Cauchy distribution.
        * 
        * @param current The current (last accepted) point
        * @param candidate The address for the candidate point to be evaluated
        * @param temperature Temperature for each dimension of the problem (these can be all different)
        */
        void     makeNew(vec &current, strc &candidate, vec &temperature);

        /**
         * @brief Get the States Since Reanneal parameter
         * 
         * @return unsigned the number of points generated since the last re-annealing event
         */
        unsigned getStatesSinceReanneal(void);

        /**
         * @brief Get the Total States parameter
         * 
         * @return unsigned the total number of points generated so far
         */
        unsigned getTotalStates(void);

       /**
        * @brief Checks if re-annealing should be performed according to the criterion of the Generator.
        * 
        * @return true if re-annealing should be performed
        * @return false otherwise
        */
        bool reAnneal(void);

        /**
         * @brief Resets the counter for the number of points generated since the last re-annealing event.
         * 
         */
        void reset(void);

       /**
        * @brief Checks if the stopping criterion for the Generator occurs.
        * 
        * @return true if the stopping criterion is satisfied
        * @return false otherwise
        */
        bool exit(void);

};
