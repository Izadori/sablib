/**
 * @file imodpoly.h
 * @brief Baseline estimation using Improved Modified Polynomial(IModPoly) method
 * @author Izadori
 * @details
 * References:
 * @li Zhao, J.; Lui, H.; McLean, D. I.; Zeng, H. Automated Autofluorescence Background Subtraction Algorithm for Biomedical Raman Spectroscopy" Applied Spectroscopy 2007, 61(11), 1225-1232.
 */

#ifndef __SABLIB_IMODPOLY_H__
#define __SABLIB_IMODPOLY_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Estimates the baseline using the Improved Modified Polynomial (IModPoly) method.
 *
 * @param y The input data points to be processed.
 * @param polyorder The order of the polynomial to fit.
 * @param k Scaling factor for the standard deviation threshold (default is 1.0).
 * @param loop The maximum number of iterations (default is 50).
 * @param eps The convergence tolerance based on the standard deviation change (default is 1.0e-3).
 * @return A vector containing the estimated baseline.
 * @throw std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double> BaselineIModPoly(
	const std::vector<double> & y, const unsigned int polyorder, const double k = 1,
	const unsigned int loop = 50, const double eps = 1.0e-3
);

}; // namespace sablib

#endif // __SABLIB_IMODPOLY_H__
