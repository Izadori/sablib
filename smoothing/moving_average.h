/**
 * @file moving_average.h
 * @brief Smoothing using simple moving average
 * @author Izadori
 */

#ifndef __SABLIB_MOVING_AVERAGE_H__
#define __SABLIB_MOVING_AVERAGE_H__

#include <cmath>
#include <vector>

#include "../misc/convolve.h"

namespace sablib {

/**
 * @brief Calculates the simple moving average of the input signal (std::vector<double> version).
 *
 * @param y The data to be averaged.
 * @param n Half-width of the moving average window (calculated using `2 * n + 1` points).
 * @return The data after applying the moving average.
 */
const std::vector<double> MovingAverage(const std::vector<double> & y, const unsigned int n);

/**
 * @brief Calculates the simple moving average of the input signal.
 *
 * @param y The data to be averaged.
 * @param n Half-width of the moving average window (calculated using `2 * n + 1` points).
 * @return The data after applying the moving average.
 */
template <typename Derived>
const typename Derived::PlainObject MovingAverage(const Eigen::MatrixBase<Derived> & y, const unsigned int n)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: y is not vector.");

	using PlainObject = typename Derived::PlainObject;

	int points = 2 * n + 1;
	PlainObject w = PlainObject::Ones(points) / points;

	PlainObject result = Convolve(y, w, ConvolveMode::Same);

	return result;
}

}; // namespace sablib

#endif // __SABLIB_MOVING_AVERAGE_H__
