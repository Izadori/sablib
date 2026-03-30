/**
 * @file psalsa.h
 * @brief Baseline estimation using Peaked Signal’s Asymmetric Least Squares Algorithm(psalsa)
 * @author Izadori
 * @details
 * References:
 * @li Oller-Moreno S.; Pardo A.; Jiménez-Soto J. M.; Samitier J.; Marco S. "Adaptive Asymmetric Least Squares baseline estimation for analytical instruments" 2014 IEEE 11th International Multi-Conference on Systems, Signals & Devices (SSD14), Barcelona, Spain, 2014, pp. 1-5.
 */

#ifndef __SABLIB_PSALSA_H__
#define __SABLIB_PSALSA_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Performs baseline estimation using Peaked Signal’s Asymmetric Least Squares Algorithm (psalsa).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param p Weight (asymmetry parameter, typically 0.001 to 0.1).
 * @param k Exponential decay of the weights.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double> BaselinePsalsa(
	std::vector<double> & y, const double lambda, const double p, const double k,
	const unsigned int s = 2, const unsigned int loop = 10, const double eps = 1e-3
);

}; // namespace sablib

#endif // __SABLIB_PSALSA_H__
