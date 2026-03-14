/**
 * @file polynomial.h
 * @brief Baseline estimation with polynomial line
 * @author Izadori
 */

#ifndef __SABLIB_POLYNOMIAL_H__
#define __SABLIB_POLYNOMIAL_H__

#include <stdexcept>
#include <vector>
#include <Eigen/Eigen>

#include "../misc/polyfit.h"

namespace sablib {

/**
 * @brief Performs baseline estimation with a linear line between two points.
 *
 * @param y The input data for baseline estimation.
 * @param index1 The index of the first point.
 * @param index2 The index of the second point.
 * @return The estimated linear baseline.
 */
std::vector<double> BaselineLinear(std::vector<double> & y, const unsigned int index1, const unsigned int index2);

/**
 * @brief Performs baseline estimation by fitting a polynomial to specified points.
 *
 * @param y The input data for baseline estimation.
 * @param polyorder The order of the polynomial to fit.
 * @param indices The indices of the points used for polynomial fitting.
 * @return The estimated polynomial baseline.
 */
std::vector<double> BaselinePolynomial(
	std::vector<double> & y, const unsigned int polyorder, const std::vector<unsigned int> & indices
);

}; // namespace sablib

#endif // __SABLIB_POLYNOMIAL_H__
