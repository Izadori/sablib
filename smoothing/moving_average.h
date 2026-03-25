/**
 * @file moving_average.h
 * @brief Smoothing using simple/weighted moving average
 * @author Izadori
 */

#ifndef __SABLIB_MOVING_AVERAGE_H__
#define __SABLIB_MOVING_AVERAGE_H__

#include <cmath>
#include <stdexcept>
#include <vector>

#include "../misc/expand.h"

namespace sablib {

/**
 * @brief Calculates the weighted moving average of the input signal (std::vector<double> version).
 *
 * @param y The data to be averaged.
 * @param w Weights.
 * @return The data after applying the moving average.
 */
const std::vector<double> WeightedMovingAverage(const std::vector<double> & y, const std::vector<double> & w);

/**
 * @brief Calculates the weighted moving average of the input signal.
 *
 * @param y The data to be averaged.
 * @param w Weights.
 * @return The data after applying the moving average.
 * @exception std::invalid_argument If the length of y or w is zero.
 */
template <typename Derived>
const typename Derived::PlainObject WeightedMovingAverage(
	const Eigen::MatrixBase<Derived> & y, const Eigen::MatrixBase<Derived> & w
)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: y and w are not vector.");

	if(y.size() == 0 || w.size() == 0) {
		throw std::invalid_argument("WeightedMovingAverage(): the length of y or w is zero.");
	}

	using PlainObject = typename Derived::PlainObject;

	int points = w.size();
	int n = points / 2;
	PlainObject yy = ExpandBoundaries(y, n);
	PlainObject result = PlainObject::Zero(y.size());

	if(y.cols() == 1) {
		for(int i = 0; i < y.size(); i++) {
			result(i) = (yy.block(i, 0, points, 1).array() * w.array()).matrix().sum();
		}
	}
	else if(y.rows() == 1) {
		for(int i = 0; i < y.size(); i++) {
			result(i) = (yy.block(0, i, 1, points).array() * w.array()).matrix().sum();
		}
	}

	return result;
}

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
 * @exception std::invalid_argument If n is zero.
 */
template <typename Derived>
const typename Derived::PlainObject MovingAverage(const Eigen::MatrixBase<Derived> & y, const unsigned int n)
{
	using PlainObject = typename Derived::PlainObject;

	if(n == 0) {
		throw std::invalid_argument("MovingAverage(): n is zero.");
	}

	int points = 2 * n + 1;
	PlainObject w = PlainObject::Ones(points) / points;

	return WeightedMovingAverage(y, w);
}

/**
 * @brief Generates a Gaussian kernel.
 *
 * @param n Half-width of the Gaussian window (total size is `2 * n + 1`).
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return A vector containing the Gaussian kernel coefficients.
 * @exception std::invalid_argument If n is zero or sigma is non-positive.
 */
const std::vector<double> GaussianKernel(const unsigned int n, const double sigma);

/**
 * @brief Performs Gaussian smoothing on the input signal.
 *
 * This is a convenience function that generates a Gaussian kernel and then
 * applies it using `WeightedMovingAverage`.
 *
 * @param y The data to be smoothed.
 * @param n Half-width of the Gaussian window (total size is `2 * n + 1`).
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The data after applying the Gaussian filter.
 */
inline const std::vector<double> GaussianFilter(const std::vector<double> & y, const unsigned int n, const double sigma)
{
	return WeightedMovingAverage(y, GaussianKernel(n, sigma));
}

}; // namespace sablib

#endif // __SABLIB_MOVING_AVERAGE_H__
