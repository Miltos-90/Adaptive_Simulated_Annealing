#include"Generator.hpp"
#include<vector>
#include<random>
#include"utils.hpp"

/*  ============================= Public Members ============================= */

// Draw candidate point from Cauchy distribution
void Generator::makeNew(vec &current, strc &candidate, vec &temperature)
{
    std::uniform_real_distribution<> U(0.0, 1.0);
    double u, dx;

    do{// Sample from the Cauchy distribution to make a new candidate point
        for(unsigned i = 0; i < noDimensions; i++) {            
            u  = U(Seed);
            dx = utils::sgn(u - 0.5) * temperature[i] * ( std::pow( 1 + 1 / temperature[i], std::abs(2 * u - 1) ) - 1 );
            candidate.variables[i] = current[i] + dx * (upperBound[i] - lowerBound[i]);

            if (isInteger(i)){
                candidate.variables[i] = std::round(candidate.variables[i]);
            }
        }
    } while( isOutOfBounds(candidate.variables) | violatesConstraints(candidate.variables));
    candidate.cost = costFunction(candidate.variables);
    reannealStates += 1;
    totalStates    += 1;

    return;
};

// Getters
unsigned Generator::getStatesSinceReanneal(void) {return reannealStates;}
unsigned Generator::getTotalStates(void) {return totalStates;}

// Set properties required for the main algorithm loop
void Generator::setProperties(bvec &intVariables, vec &lbound, vec &ubound, fcn &func, bfcn &constrFunction, std::minstd_rand &rEngine)
{
    lowerBound          = lbound;
    upperBound          = ubound;
    integerVars         = intVariables;
    costFunction        = func;
    constraintFunction  = constrFunction;
    noDimensions        = lbound.size();
    Seed                = rEngine;
}

// Checks if a generated point violates any constraints
bool Generator::violatesConstraints(vec &point)
{
    std::vector<bool> contraints = constraintFunction(point);
    return !std::all_of(std::begin(contraints), std::end(contraints),  [](bool i) {return i;});;
}

// Reset counter
void Generator::reset(void) { reannealStates = 0; return; }

// Check if it's time for reAnnealing according to the generated points criterion
bool Generator::reAnneal(void){return reannealStates == frequencyModulus;}

// Check if it's time to exit according to the total generated states
bool Generator::exit(void){return totalStates <= maxStates;}

/*  ============================ Private Members ============================= */

// Check if drawn point lies within the search space
bool Generator::isOutOfBounds(vec &point)
{
    bool outOfBounds = false; // Assume point is feasible

    // Loop over all dimensions
    for (unsigned i = 0; i < point.size(); i ++)
    {
        // Check if point dimension is outside of bounds
        if ( (point[i] < lowerBound[i]) | (point[i] > upperBound[i]) )
        {
            outOfBounds = true;
            break; // Exit if unfeasible dimension is found
        }
    }

    return outOfBounds;
}

// Check if an input dimension of the problem corresponds to an integer variable
bool Generator::isInteger(unsigned dimension){ return integerVars[dimension];}