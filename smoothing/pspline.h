/**
 * @file pspline.h
 * @brief Smoothing using penalized Spline(P-Spline)
 * @author Izadori
 * @details
 * References:
 * @li Eilers, P. H. C.; Marx, B. D. "Flexible Smoothing with B-Splines and Penalties" Stastical Science, 1996, 11(2), 89-121.
 */

#ifndef __SABLIB_PSPLINE_H__
#define __SABLIB_PSPLINE_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Smoothes the input data using P-Splines (Penalized B-Splines).
 *
 * This function applies a B-Spline basis with a difference penalty on the coefficients
 * to achieve smoothing of the input data.
 *
 * @param y The input data points to be smoothed.
 * @param knots_num The number of internal knots for the B-Spline basis.
 * @param degree The degree of the B-Spline basis (default is 3, cubic spline).
 * @param s The order of the difference penalty (default is 2).
 * @param lambda The smoothing parameter (default is 1.0). Larger values result in smoother curves.
 * @return A vector containing the smoothed data points.
 * @throw std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double> PSpline(
	const std::vector<double> & y, const unsigned int knots_num,
	const unsigned int degree = 3, const unsigned int s = 2, const double lambda = 1.0
);

}; // namespace sablib

#endif // __SABLIB_PSPLINE_H__
