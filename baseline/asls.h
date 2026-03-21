/**
 * @file asls.h
 * @brief Baseline estimation using Asymmetric Least Squares Smoothing(AsLS)
 * @author Izadori
 * @details
 * References:
 * @li Eilers, P. A Perfect Smoother. Analytical Chemistry, 2003, 75(14), 3631-3636.
 * @li Eilers, P. H. C.; Boelens, H. F. M. Baseline Correction with Asymmetric Least Squares Smoothing Leiden University Medical Centre Report 2005, 1(1), 1-5.
 */

#ifndef __SABLIB_ASLS_H__
#define __SABLIB_ASLS_H__

#include <stdexcept>

#include "../smoothing/whittaker.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using Asymmetric Least Squares Smoothing (AsLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param p Weight (asymmetry parameter, typically 0.001 to 0.1).
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double> BaselineAsLS(
	std::vector<double> & y, const double lambda, const double p, const unsigned int s = 2,
	const unsigned int loop = 10, const double eps = 1e-3
);

}; // namespace sablib

#endif // __SABLIB_ASLS_H__
