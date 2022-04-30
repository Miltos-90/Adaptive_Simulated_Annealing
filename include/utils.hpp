/**
 * @file utils.hpp
 * @author miltos-90
 * @brief tils namespace. Implements some genetic templates used throughout.
 * @version 1.0
 * @date 2022-04-17
 * 
 */

#pragma once
#include <vector>
#include <functional>
#include <limits>
#include<cmath>

namespace utils
{
    /**
     * @brief Type-safe signum function
     * 
     * @tparam T datatype of the input
     * @param val input value
     * @return int the sign of the input 
     */
    template <typename T> 
    int sgn(T val) { return (T(0) < val) - (val < T(0)); }

    // 

    /**
     * @brief Struct to hold decision variables and corresponding cost function
     * 
     * @tparam T the datatype
     */
    template <typename T> 
    struct pointStruct{
        // Fields
        std::vector<T> variables; 
        T cost = std::numeric_limits<T>::max();
        // Constructor
        pointStruct(size_t vSize): variables(vSize, 0.0){};
        };
}