/**
 * @file ASA.hpp
 * @author miltos-90
 * @brief Header file for the ASA (Adaptive Simulated Annealing) algorithm.
 * 
 * This is a metaheuristic that can handle singe objective mixed-integer contrained optimisation problems of the form:
 * min(x)  f(x)
 * s.t.:   G(x)
 *          LowerBound_i <= x_i <= UpperBound_i
 *          (optionally) x_i in Z,
 *          for all i {1,..,Nx}
 * where:
 *          x:    vector<double> with <noDimensions> elements
 *          f(x): double
 *          G(x): vector<bool> with <NumConstraints> elements (=True if a constraint is satisfied / =False is constraint is violated)
 * The algorithm is described in:
 * [1] Ingber, L. (1989). Very fast simulated re-annealing. Mathematical and computer modelling, 12(8), 967-973.
 * [2] Ingber, L. (1993). Adaptive simulated annealing (ASA). Global optimization C-code, Caltech Alumni Association, Pasadena, CA.
 * [3] Ingber, L. (1993). Simulated annealing: Practice versus theory. Mathematical and computer modelling, 18(11), 29-57.
 * [4] Ingber, L. (2000). Adaptive simulated annealing (ASA): Lessons learned. arXiv preprint cs/0001018.
 * 
 * @version 1.0
 * @date 2022-04-16
 * 
 */

#pragma once
#include<vector>
#include<functional>
#include<random>
#include<thread>
#include<stdexcept>
#include<iostream>
#include<stdlib.h>
#include"utils.hpp"
#include"Generator.hpp"
#include"Assessor.hpp"
#include"TemperatureScheduler.hpp"
#include"ThreadPool.hpp"
#include"Assessor.hpp"

class ASA{
    using vec  = std::vector<double>;
    using bvec = std::vector<bool>;
    using strc = utils::pointStruct<double>;                // Struct containing a point in the N-dimensional input space and it associated cost value.
    using fcn  = std::function<double (vec &x)>;            // Type of the cost (objective) function
    using bfcn = std::function<std::vector<bool> (vec &x)>; // Type of the contraint function

    private:
        
        // Fields
        unsigned  costRepeats = 0;          // Counter for the no. times |f_best - f_previousBest| <= costEps so far
        double    costEps;                  // Max. cost difference assumed to be no improvement
        unsigned  maximumCostRepeat;        // Max no. times |f_best - f_previousBest| <= costEps (stopping criterion)
        double    acceptedToGeneratedRatio; // Ratio of accepted to generated states (stopping criterion)
        Generator PointGen;                 // Object to generate candidate points
        Assessor  PointAss;                 // Object to assess candidate points
        TemperatureScheduler Scheduler;     // Object to regulate the annealing temperatures
        
        // Members
       
       /**
        * @brief Implements the ASA algorithm.
        * 
        * @param Sa An object of this class (to support multi-threading)
        * @param lowerBound Vector containing the lower bound of the search space for each dimension
        * @param upperBound Vector containing the upper bound of the search space for each dimension
        * @param intVariables Vector indicating which decision variables are continous/integers
        * @param costFunction Objective function (double Vector)->double that will be minimised
        * @param constraintFunction Contstraints(double Vector) -> Bool vector indicating if each constraint is satisfied (=1) or violated (=0)
        * @return strc struct containing the best obtained point and the value of its cost function
        */
        static strc main(ASA Sa, vec &lowerBound, vec &upperBound, bvec &intVariables,
                         fcn &costFunction, bfcn &constraintFunction);

       /**
        * @brief updates the best point found so far, and the cost of the previous best point found so far (needed for stopping criteria).
        * 
        * @param currentPnt Struct containing the point generated at the current iteration
        * @param bestPnt Struct containing the best point found so far by the algorithm
        * @param previousBestCost The cost (objective function value) of the previous best point found so far
        * @return double The updated value of the previousBestCost
        */
        static double updateBest(strc &currentPnt, strc &bestPnt, double previousBestCost);

       /**
        * @brief Implements the re-annealing of all temperatures and counters related to the decision variables and the cost
        * 
        * @param currentPnt Struct containing the point generated at the current iteration
        * @param bestCost The cost (objective function value) of the best point found so far
        * @param previousBestCost The cost of the previous best point found so far
        */
        void reAnneal(strc &currentPnt, double bestCost, double previousBestCost);

        /**
         * @brief Checks the stopping critera at each iteration of the main algorithm 
         * 
         * @return true if the algorithm should exit
         * @return false for the next iteration
         */
        bool stoppingCriteria(void);

       /**
        * @brief Checks if the number of accepted states / number of generatred states is lower than a threshold value set by the user (stopping criterion).
        * 
        * @return true it the criterion is satisfied
        * @return false otherwise
        */
        bool reAnnealRatioCriterion(void);

        /**
         * @brief Implements a basic error check on the inputs. Throws excpetion if an error is detection
         * 
         * @param lowerBound Vector containing the lower bound of the search space for each dimension
         * @param upperBound Vector containing the upper bound of the search space for each dimension
         * @param intVariables Vector indicating which decision variables are continous/integers
         */
        void inputCheck(vec &lowerBound, vec &upperBound, bvec &intVariables);
    
    public:
    
       /**
        * @brief Constructs a new ASA object with default parameters
        * 
        * @param TemperatureRatioScale This scale is a guide to the expected cost temperature of convergence within a small range of the global minimum.
        * @param TemperatureAnnealScale This scale is a guide to achieve the expected cost temperature sought by Temperature Ratio Scale
        * @param costScaleRatio This is the ratio of cost to parameter temperature annealing scales
        * @param acceptanceToGenerationRatio The least ratio of accepted to generated states. If this value is encountered, then the usual tests, including 
        * possible reannealing, are initiated evenifthe timing does not coincide with Acceptance Frequency Modulus or Generated Frequency Modulus
        * @param costPrecision sets the precision required of the cost function if exiting because of reaching Maximum Cost Repeatitions
        * @param generationFrequencyModulus The frequency of testing for periodic testing and reannealing, dependent on the number of generated states
        * @param acceptanceFrequencyModulus The frequency of testing for periodic testing and reannealing, dependent on the number of accepted states
        * @param maximumCostRepetitions The maximum number of times that the cost function repeats itself, within limits set by Cost Precision, before quitting.
        * @param maxGeneratedStates The maximum number of states generated before quitting.
        * @param maxAcceptedStates The maximum number of states accepted before quitting.
        * @param numberCostSamples Number of samples used to estimate the average initial cost temperature.
        */
        ASA(double   TemperatureRatioScale       = 1e-5, 
            double   TemperatureAnnealScale      = 1e2, 
            double   costScaleRatio              = 1.0,
            double   acceptanceToGenerationRatio = 1e-2,
            double   costPrecision               = 1e-18,
            unsigned generationFrequencyModulus  = 1e4, 
            unsigned acceptanceFrequencyModulus  = 1e2, 
            unsigned maximumCostRepetitions      = 10,
            unsigned maxGeneratedStates          = 99999,
            unsigned maxAcceptedStates           = 5e4,
            unsigned numberCostSamples           = 30
            ):
                costEps(costPrecision),
                maximumCostRepeat(maximumCostRepetitions),
                acceptedToGeneratedRatio(acceptanceToGenerationRatio),  
                PointGen(Generator(generationFrequencyModulus, maxGeneratedStates)),
                PointAss(Assessor(acceptanceFrequencyModulus,  maxAcceptedStates)),
                Scheduler(TemperatureScheduler(TemperatureRatioScale, TemperatureAnnealScale, costScaleRatio, numberCostSamples))
                {};
        
       /**
        * @brief The function called by the user to setup and solve the given optimisation (minimisation) problem.
        * 
        * @param costFunction The objective function to be minimised
        * @param lowerBound The lower bounds of the decision variables
        * @param upperBound The upper bounds of the decision variables
        * @param integerVariables Indicator for continuous/integer devision variables
        * @param constrFunction Function containing the constraints to be satisfied
        * @param numThreads Number of threads to be used for the multi-starts
        * @param numMultistarts Number of times to re-start the algorithm (with different initial conditions)
        *
        * @return vec the optimum point
        */
        vec minimize(fcn costFunction, vec &lowerBound, vec &upperBound, bvec &integerVariables,
                      unsigned numThreads = 0, unsigned numMultistarts = 1,
                      bfcn constrFunction = [](vec &) -> bvec{ return bvec(1, true); });

};