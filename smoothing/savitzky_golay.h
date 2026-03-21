/**
 * @file savitzky_golay.h
 * @brief Smoothing using Savitzky-Golay filter
 * @author Izadori
 * @details
 * References:
 * @li Savitzky, A.; Golay, M. J. Smoothing And differentiation of Data by Simplified Least Squares Procedures. Anal. Chem. 1964, 36(8), 1627-1639.
 */

#ifndef __SABLIB_SAVITZKY_GOLAY_H__
#define __SABLIB_SAVITZKY_GOLAY_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Calculates the coefficients for a Savitzky-Golay filter.
 *
 * @param n The half-width of the filter window (total window size is `2 * n + 1`).
 * @param polyorder The order of the polynomial used for fitting.
 * @param derive The order of the derivative to compute (0 for smoothing only).
 * @param delta The spacing of the samples.
 * @return The calculated filter coefficients.
 * @throw std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double> SavitzkyGolayCoefficients(
	const unsigned int n, const unsigned int polyorder, const unsigned derive = 0, const double delta = 1
);

/**
 * @brief Performs smoothing (and differentiation) using a Savitzky-Golay filter.
 *
 * @param y The input data to be filtered.
 * @param n The half-width of the filter window (total window size is `2 * n + 1`).
 * @param polyorder The order of the polynomial used for fitting.
 * @param derive The order of the derivative to compute (0 for smoothing only).
 * @param delta The spacing of the samples.
 * @return The filtered (smoothed or differentiated) data.
 * @throw std::invalid_argument If the length of y is zero.
 */
const std::vector<double> SavitzkyGolay(
	const std::vector<double> & y,
	const unsigned int n, const unsigned int polyorder, const unsigned derive = 0, const double delta = 1
);

}; // namespace sablib

#endif // __SABLIB_SAVITZKY_GOLAY_H__
