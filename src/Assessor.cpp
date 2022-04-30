#include"Assessor.hpp"
#include<vector>
#include<random>
#include"utils.hpp"

/*  ============================= Public Members ============================= */

// Accepts a point
void Assessor::accept(strc &currentPoint, strc &candidatePoint, double costTemperature)
{
    // Compute Boltzmann test
    std::uniform_real_distribution<> U(0.0, 1.0);
    bool bTest = std::exp( - (candidatePoint.cost - currentPoint.cost ) / costTemperature ) > U(Seed);
    
    // If yes, swap it and increment counters
    if (bTest){
        reannealStates += 1;
        totalStates    += 1;
        currentPoint   = candidatePoint;
    }
    
    return;
}

// Initialise properties
void Assessor::setProperties(std::minstd_rand &rEngine)
{
    Seed = rEngine;
    return;
}

// Reset counter for every re-anneal event
void Assessor::reset(void){reannealStates = 0; return;}

// Getters
unsigned Assessor::getStatesSinceReanneal(void) {return reannealStates;}
unsigned Assessor::getTotalStates(void) {return totalStates;}

// Check if it's time for reAnnealing according to the accepted points criterion
bool Assessor::reAnneal(void){return reannealStates == frequencyModulus;}

// Check if it's time to exit according to the generated states criterion 
bool Assessor::exit(void){return totalStates < maxStates;}
