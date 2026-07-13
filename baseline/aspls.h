/**
 * @file aspls.h
 * @brief Baseline estimation using adaptive smoothness Penalized Least Squares(asPLS)
 * @author Izadori
 * @details
 * References:
 * @li Zhang, F.; Tang, X.-J.; Tong, A., et al. "Baseline correction for infrared spectra usingadaptive smoothness parameter penalized leastsquares method" Spectroscopy Letters 2020, 53(3), 1-12.
 */

#ifndef __SABLIB_ASPLS_H__
#define __SABLIB_ASPLS_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using adaptive smoothness Penalized Least Squares(asPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param k Asymmetric coefficient.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselineAsPLS(
	const std::vector<double> & y, const double lambda, const unsigned int k = 0.5,
	const unsigned int s = 2, const unsigned int loop = 50, const double eps = 1e-3
);

}; // namespace sablib

#endif // __SABLIB_ASPLS_H__
