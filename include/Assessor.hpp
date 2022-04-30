/**
 * @file Assessor.hpp
 * @author miltos-90
 * @brief Header file for the Assessor class. 
 * 
 * Header file for the Assessor class. This class implements the so-called Boltzmann test
 * to accept newly-generated points during an iteration of the ASA algorithm.
 * Used in the ASA class.
 * 
 * @version 1.0
 * @date 2022-04-17
 * 
 */

#pragma once
#include<vector>
#include<random>
#include"utils.hpp"

class Assessor
{
    using vec  = std::vector<double>;
    using strc = utils::pointStruct<double>;

    private:
        // Fields
        unsigned reannealStates = 0; // Counter for the number of points accepted since the last re-anneling
        unsigned totalStates    = 0; // Counter for the total number of accepted points so far
        unsigned frequencyModulus;   // Frequency for re-annealing
        unsigned maxStates;          // Max points to be assessed (stopping criteria)
        std::minstd_rand Seed;       // Seeder for random distributions
    
    public:
        
        /**
         * @brief Construct a new Assessor object
         * 
         */
        Assessor(){};
        Assessor(unsigned modulus, unsigned states):
                    frequencyModulus(modulus), maxStates(states){};

        // Members

        /**
        * @brief Implements the Boltzmann test to accept a newly-generated point (from the Generator).
        * 
        * @param currentPnt The current (last accepted) point
        * @param candidatePnt The point created by the Generator class
        * @param costTemperature The temperature related to the cost function, computed from the Temperature Scheduler
        */
        void accept(strc &currentPnt, strc &candidatePnt, double costTemperature);

        /**
         * @brief Sets the non-initialised properties
         * 
         * @param rEngine Seeder for the random distributions
         */
        void setProperties(std::minstd_rand &rEngine);

        /**
         * @brief Get the States Since Reanneal parameter
         * 
         * @return unsigned the number of points accepted since the last re-annealing event
         */
        unsigned getStatesSinceReanneal(void);

        /**
         * @brief Returns the total states parameter
         * 
         * @return unsigned the total number of points accepted so far
         */
        unsigned getTotalStates(void);

       /**
        * @brief  Checks if re-annealing should occur according to the Assessor criterion
        * 
        * @return true if the number of generated points since the previous re-annealing event exceeds the threshold given by the user
        * @return false otherwise
        */
        bool reAnneal(void);

        /**
         * @brief Resets the counter for the number of points accepted since the lst re-annealing event
         * 
         */
        void reset(void);

        /**
         * @brief Checks if the stopping criterion for the Assessor is satisfied or not.
         * 
         * @return true if the criterion is satisfied
         * @return false otherwise
         */
        bool     exit(void);
};
