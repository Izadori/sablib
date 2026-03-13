/**
 * @file arpls.h
 * @brief Baseline estimation using asymmetrically reweighted Penalized Least Squares(arPLS)
 * @author Izadori
 * @details
 * References:
 * @li Baek, S.-J.; Park, A.; Ahn, Y.-J.; Choo, J. Baseline Correction Using Asymmetrically Reweighted Penalized Least Squares Smoothing. Analyst 2015, 140(1), 250–257.
 */

#ifndef __SABLIB_ARPLS_H__
#define __SABLIB_ARPLS_H__

#include "../smoothing/whittaker.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using asymmetrically reweighted Penalized Least Squares(arPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline.
 */
const std::vector<double> BaselineArPLS(
	std::vector<double> & y, const double lambda, const unsigned int s = 2,
	const unsigned int loop = 50, const double eps = 1e-3
);

}; // namespace sablib

#endif // __SABLIB_ARPLS_H__
