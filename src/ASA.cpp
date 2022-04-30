#include"ASA.hpp"
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

/*  ============================= Public Members ============================= */

// Sets up and solves the minimisation problem
ASA::vec ASA::minimize(fcn objFunction, vec &boundLow, vec &boundHigh, bvec &intVars, 
                        unsigned numThreads, unsigned numMultistarts, bfcn constrFunction)
{
    // Input check
    try{ 
        inputCheck(boundLow, boundHigh, intVars); 
    } catch(std::invalid_argument& e) {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    // Setup some values
    unsigned maxThreads = std::thread::hardware_concurrency();
    numThreads          = (numThreads <= maxThreads) ? numThreads : maxThreads; // Cap threads
    numMultistarts      = (numMultistarts <= 1) ? 1 : numMultistarts;           // Check multistarts input
    std::vector<strc>   results(numMultistarts, 0.0);                           // Vector to hold multistart results
    strc                bestPoint(boundLow.size());                             // Best point found

    switch (numThreads) {
        case 0: // Serial computation
            for (unsigned i = 0; i < numMultistarts; i++) {
                results.push_back(main(*this, boundLow, boundHigh, intVars, objFunction, constrFunction));
            }
        break;
        default: // Multithreaded
            thread_pool pool(numThreads);
            for (unsigned i = 0; i < numMultistarts; i++) {
                pool.push_task(
                    [&results, i, this, &boundLow, &boundHigh, &intVars, &objFunction, &constrFunction] {
                        results[i] = main(*this, boundLow, boundHigh, intVars, objFunction, constrFunction);
                        });
            }
            pool.wait_for_tasks();
            break;
    }
    
    // Search the best solution among all runs
    for (strc result:results){
        if (result.cost < bestPoint.cost) { bestPoint = result; }
    }

   return bestPoint.variables;
}

/*  ============================ Private Members ============================= */

// Main implementation of the adaptive simulated annealing algorithm
ASA::strc ASA::main(ASA Sa, vec &lowerBound, vec &upperBound, bvec &intVariables, 
                    fcn &costFunction, bfcn &constraintFunction)
{
    // Fields (passed mostly by reference everywhere) that all functions operate on
    double previousBestCost = std::numeric_limits<double>::max();
    std::minstd_rand thread_local Seed(std::random_device{}());
    size_t noDimensions = lowerBound.size();
    strc   bestPoint(noDimensions);
    strc   candidatePoint(noDimensions);
    strc   currentPoint(noDimensions);

    // Initialise properties of the various objects
    Sa.Scheduler.setProperties(intVariables, lowerBound, upperBound, costFunction, Seed);
    Sa.PointGen.setProperties(intVariables, lowerBound, upperBound, costFunction, constraintFunction, Seed);
    Sa.PointAss.setProperties(Seed);
    
    // Set current point to lower bound
    currentPoint.variables  = lowerBound; 
    currentPoint.cost       = costFunction(currentPoint.variables);

    // Main loop
    double *costTemperature = Sa.Scheduler.getCostTemperature();
    vec    *temperatures    = Sa.Scheduler.getTemperatures();
    do {
        Sa.PointGen.makeNew(currentPoint.variables, candidatePoint, *temperatures); 
        Sa.PointAss.accept(currentPoint, candidatePoint, *costTemperature);
        Sa.Scheduler.step();
        Sa.reAnneal(currentPoint, bestPoint.cost, previousBestCost);
        previousBestCost = Sa.updateBest(currentPoint, bestPoint, previousBestCost);

    } while(Sa.stoppingCriteria());

    return bestPoint;
}

// Implements a basic input check
void ASA::inputCheck(vec &lowerBound, vec &upperBound, bvec &intVariables){

    // Check if sizes of vectors match
    bool sizeEqual = (upperBound.size() == lowerBound.size()) & (upperBound.size() == intVariables.size());
    
    // Check if lower bound is always lower than upper bound
    if (sizeEqual){
        for (unsigned dim = 0; dim < upperBound.size(); dim++) {
            if (upperBound[dim] <= lowerBound[dim]) {
                throw std::invalid_argument( "Detected upper bound <= lower bound. exiting" );
            }
        }
    } else {
        throw std::invalid_argument( "Unequal vector sizes detected. exiting" );
    }

    return;
}

// Keep track of the best point
double ASA::updateBest(strc &currentPnt, strc &bestPnt, double previousBestCost)
{
    if (currentPnt.cost < bestPnt.cost) {
        previousBestCost = bestPnt.cost;
        bestPnt = currentPnt;
    }
    return previousBestCost;
}

// Check stopping criteria
bool ASA::stoppingCriteria(void)
{
    return (costRepeats < maximumCostRepeat) & PointGen.exit() & PointAss.exit();
}

// Compute re-annealing criterion according to the ratio of accepted / generated points
bool ASA::reAnnealRatioCriterion(void)
{
    unsigned genPnts   = PointGen.getStatesSinceReanneal();
    unsigned accPnts   = PointAss.getStatesSinceReanneal();
    bool     condition = (accPnts > 0) & (genPnts > 0) & 
                         ( (double) accPnts / (double) genPnts <= acceptedToGeneratedRatio);   
    return condition;
}

// Update everything for the next iteration
void ASA::reAnneal(strc &currentPoint, double bestCost, double previousBestCost)
{  
    // Check if it's time for re-annealing
    if (reAnnealRatioCriterion() | PointGen.reAnneal() | PointAss.reAnneal()){
        
        Scheduler.reAnneal(bestCost, currentPoint);

        // Update cost repeat counter
        if (std::abs(previousBestCost - bestCost) < costEps) {costRepeats += 1;}

        // Reset counters
        PointAss.reset();
        PointGen.reset();
    }
    
    return;
}
