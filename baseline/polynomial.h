/**
 * @file polynomial.h
 * @brief Baseline estimation with polynomial line
 * @author Izadori
 */

#ifndef __SABLIB_POLYNOMIAL_H__
#define __SABLIB_POLYNOMIAL_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs baseline estimation with a linear line between two points.
 *
 * @param y The input data for baseline estimation.
 * @param index1 The index of the first point.
 * @param index2 The index of the second point.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselineLinear(const std::vector<double> & y, const unsigned int index1, const unsigned int index2);

/**
 * @brief Performs baseline estimation by fitting a polynomial to specified points.
 *
 * @param y The input data for baseline estimation.
 * @param polyorder The order of the polynomial to fit.
 * @param indices The indices of the points used for polynomial fitting.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselinePolynomial(
	const std::vector<double> & y, const unsigned int polyorder, const std::vector<unsigned int> & indices
);

}; // namespace sablib

#endif // __SABLIB_POLYNOMIAL_H__
