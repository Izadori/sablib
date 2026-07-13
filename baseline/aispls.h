/**
 * @file aispls.h
 * @brief Baseline estimation using adaptive iterative smoothing parameter Penalized Least Squares(aisPLS)
 * @author Izadori
 * @details
 * References:
 * @li He, J.; Bai, Y.; Jv, Z.; Chen, Z.; Wang, B. "Adaptive Multi-Order Penalty and Dual-Driven Weighting: aisPLS Algorithm for Raman Baseline Correction with Weak Peak Preservation" Molecules 2026, 31(8), 1243–1258.
 */

#ifndef __SABLIB_AISPLS_H__
#define __SABLIB_AISPLS_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using adaptive iterative smoothing parameter Penalized Least Squares(aisPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param r Adaptation rate of smoothing parameter.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselineAisPLS(
	const std::vector<double> & y, const double lambda, const double r = 2,
	const unsigned int s = 2, const unsigned int loop = 50, const double eps = 1e-3
);

}; // namespace sablib

#endif // __SABLIB_AISPLS_H__
