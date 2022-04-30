#include"TemperatureScheduler.hpp"
#include"Derivatives.hpp"
#include"utils.hpp"
#include<algorithm>
#include<vector>
#include<cmath>
#include<random>

/*  ============================= Public Members ============================= */

// Increment all counters and temperatures after one iteration of the main algorithm
void TemperatureScheduler::step(void)
{
    // Update cost-related counters and temperatures
    costIterationCounter += 1;
    currentCostTemperature = initialCostTemperature * std::exp(- costParameterRatio * std::pow(costIterationCounter, invDimensions) );
    
    // Update decision variable temperatures and counters
    for (unsigned dim = 0; dim < noDimensions; dim++){ 
        iterationCounters[dim]   = iterationCounters[dim] + 1;
        currentTemperatures[dim] = initialTemperatures[dim] * std::exp(-temperatureParameterRatio * std::pow(iterationCounters[dim], invDimensions) );
    }
    return;
}

// Reanneal: Collect all reannealing-related functions in a single call
void TemperatureScheduler::reAnneal(double bestCost, strc &currentPoint)
{
    reAnnealInitialCostTemperature(bestCost, currentPoint.cost);
    reAnnealCostIterationCounter(bestCost, currentPoint.cost);
    reAnnealTemperatures(currentPoint.variables);
    return;
}

// Setter for various properties needed for the main loop of simulated annealing
void TemperatureScheduler::setProperties(bvec &intVariables, vec &lowerBound, vec &upperBound, fcn &func, std::minstd_rand &rEngine)
{
    diff                      = Derivatives(lowerBound, upperBound, func);
    noDimensions              = lowerBound.size();
    temperatureParameterRatio = m * std::exp(-n / (double) noDimensions);
    costParameterRatio        = temperatureParameterRatio * costParameterScaleRatio;
    initialTemperatures       = vec(noDimensions, 1.0);
    iterationCounters         = vec(noDimensions, 0.0);
    sensitivities             = vec(noDimensions, 0.0);
    currentTemperatures       = initialTemperatures;
    invDimensions             = 1.0 / (double) noDimensions;
    costFunction              = func;
    integerVars               = intVariables;
    getInitialTemperature(numberCostSamples, lowerBound, upperBound, rEngine);  
    return;
}

// Various getters
double* TemperatureScheduler::getCostTemperature(){ return &currentCostTemperature; }
std::vector<double>* TemperatureScheduler::getTemperatures(){ return &currentTemperatures; }

/*  ============================= Private Members ============================= */

// Check if an input dimension of the problem corresponds to an integer variable
bool TemperatureScheduler::isInteger(unsigned dimension){ return integerVars[dimension];}

// Estimate initial cost temperature from a set of random function calls
void TemperatureScheduler::getInitialTemperature(unsigned noRepeats, vec &lowerBound, vec &upperBound, std::minstd_rand &rEngine)
{
    unsigned i, j;
    vec point(noDimensions, 0.0);
    std::uniform_real_distribution<> U(0.0, 1.0);

    for (i = 0; i < noRepeats; i++) {

        // Draw a random point within the admissible range
        for (j = 0; j < noDimensions; j++) {
            point[j] = lowerBound[j] + U(rEngine) * (upperBound[j] - lowerBound[j]);

            if (isInteger(j)){
                point[j] = std::round(point[j]);
            }
        }

        // Compute cost
        initialCostTemperature += std::abs(costFunction(point)) / (double) noRepeats; 
    }
    // Set current cost temperature
    currentCostTemperature = initialCostTemperature;

    return;
}

// Reanneal initial cost temperature
void TemperatureScheduler::reAnnealInitialCostTemperature(double bestCost, double currentCost)
{
    double maxCost         = std::max(std::abs(bestCost), std::abs(currentCost));
    double deltaCost       = std::abs(bestCost - currentCost);
    double eps             = std::numeric_limits<double>::min();
    double maxF            = std::max(eps, std::max(maxCost, deltaCost));
    initialCostTemperature = std::min(initialCostTemperature, maxF);

    return;
}

// Reanneal cost iteration counters
void TemperatureScheduler::reAnnealCostIterationCounter(double bestCost, double currentCost)
{
    double deltaCost     = std::abs(bestCost - currentCost);
    double tmp           = std::max(deltaCost, initialCostTemperature);
    double maxF          = std::max(std::numeric_limits<double>::min(), tmp);
    double reAnnCostT    = std::min(initialCostTemperature, maxF);
    double tRatio        = initialCostTemperature / reAnnCostT;
    costIterationCounter = (unsigned) std::pow(std::log10(tRatio) / costParameterRatio, (double) noDimensions);
    
    return;
}

// Reanneal decision variable temperatures
void TemperatureScheduler::reAnnealTemperatures(vec &point)
{
    double derivative, sMax = std::numeric_limits<double>::min();
    double newTemperature, logTemperatureRatio;
    
    // Compute sensitivities (partial derivatives) and assign the maximum 
    for (unsigned k = 0; k < noDimensions; k++) {
        if (!isInteger(k)) {
            derivative = std::abs(diff.compute(point, k));
            if (derivative > sMax) {sMax = derivative;}
            sensitivities[k] = derivative;
        }
    }

    // Compute new iteration counters
    for (unsigned k = 0; k < noDimensions; k++) {
        if (!isInteger(k)) {    
            newTemperature       = std::min(initialTemperatures[k], currentTemperatures[k] * sMax / sensitivities[k]);
            logTemperatureRatio  = std::log10(initialTemperatures[k] / newTemperature);
            iterationCounters[k] = std::round(std::pow(logTemperatureRatio / temperatureParameterRatio, (double) noDimensions));
        }
    }

    return;
}
