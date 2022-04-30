#include <iostream>
#include <cmath>
#include "ASA.hpp"

/* =================================================================================
    Minimisation example of the Rosenbrock function constrained with a cubic and a line
    (https://en.wikipedia.org/wiki/Test_functions_for_optimization)
    using the Adaptive Simulated Annealing Algorithm.
 ===================================================================================*/

// Objective function
double Rosenfunction(std::vector<double> &x){
    /*
    This function computes what is known as Rosenbrock's function.  It is 
    a function of two input variables and has a global minimum at (1,1).
    So when we use this function to test out the optimization algorithms
    we will see that the minimum found is indeed at the point (1,1). 
    */
    return std::pow( (1 - x[0]), 2.0 ) + 100 * std::pow((x[1] - x[0] * x[0]), 2.0);
}

// Constraints function
std::vector<bool> constraints(std::vector<double> &x){ 
    /*
    This is an example constraints function that can be optionally 
    included in the optimiser. If included, it should return an 
    std::vector<bool> (numberOfContraints, true/false) indicating
    if each constraint is satisfied.
    */
    std::vector<bool> constraintVec(2, true);
    constraintVec[0] = std::pow(x[0] - 1.0, 3.0) - x[1] + 1.0 <= 0;
    constraintVec[1] = x[0] + x[1] - 2.0 <= 0;

    return constraintVec;
}

int main(void)
{
    ASA Sa/*(double   TemperatureRatioScale      = 1e-5, 
            double   TemperatureAnnealScale      = 1e2, 
            double   costScaleRatio              = 1.0,
            double   acceptanceToGenerationRatio = 1e-2,
            double   costPrecision               = 1e-18,
            unsigned generationFrequencyModulus  = 1e4, 
            unsigned acceptanceFrequencyModulus  = 1e2, 
            unsigned maximumCostRepetitions      = 10,
            unsigned maxGeneratedStates          = 99999,
            unsigned maxAcceptedStates           = 5e4,
            unsigned numberCostSamples           = 30)
          */
        ;// Initialise ASA object (w/ default parameters here)
    
    std::vector<double> lowBound = {-1.5, -0.5};    // Set lower bound
    std::vector<double> upBound  = {1.5, 2.5};      // Set upper bound
    std::vector<bool>   intVars  = {false, false};  // Set a boolean vector indicating integer variables
    unsigned threads     = 5;                       // Number of threads to be used
    unsigned multistarts = 5;                       // Number of multistarts to be performed

    // Minimise the problem
    std::vector<double> xOpt = 
        Sa.minimize(Rosenfunction,  // Objective
                    lowBound,       // Lower bound
                    upBound,        // Upper bound
                    intVars,        // Integer constraints
                    threads,        // Number of threads
                    multistarts,    // Multistarts
                    constraints     // linear & non-linear constraints (comment out if not applicable)
                    );
    
    // Print result
    std::cout << "Optimal point: ";
    for(double v: xOpt) { std::cout << v << ", "; }
    std::cout << "Optimal cost: " << Rosenfunction(xOpt) << std::endl;
}