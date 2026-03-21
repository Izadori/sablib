/**
 * @file airpls.h
 * @brief Baseline estimation using adaptive iteratively reweighted Penalized Least Squares(airPLS)
 * @author Izadori
 * @details
 * References:
 * @li Zhang, Z.-M.; Chena, Shan; Liang, Y.-Z. Baseline correction using adaptive iteratively reweighted penalized least squares Analyst 2010, 135(5), 1138-1146.
 * @li <a href="https://github.com/zmzhang/airPLS">Zhang's GitHub</a>
 */

#ifndef __SABLIB_AIRPLS_H__
#define __SABLIB_AIRPLS_H__

#include <stdexcept>

#include "../smoothing/whittaker.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using adaptive iteratively reweighted Penalized Least Squares(airPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double> BaselineAirPLS(
	std::vector<double> & y, const double lambda, const unsigned int s = 2,
	const unsigned int loop = 50, const double eps = 1e-3
);

}; // namespace sablib

#endif // __SABLIB_AIRPLS_H__
