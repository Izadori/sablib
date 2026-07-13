/**
 * @file functions.h
 * @brief Mathematical functions for baseline estimation
 * @author Izadori
 */

#ifndef __SABLIB_FUNCTIONS_H__
#define __SABLIB_FUNCTIONS_H__

#define _USE_MATH_DEFINES
#include <cmath>

namespace sablib {

/**
 * @brief Calculates the logistic (sigmoid) function.
 *
 * Computes f(x) = 1 / (1 + e^(-x)).
 *
 * @param x The input value.
 * @return The logistic function evaluated at x.
 */
const double logistic(const double x);

} // namespace sablib

#endif // __SABLIB_FUNCTIONS_H__
