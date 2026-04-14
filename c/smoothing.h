/**
 * @file smoothing.h
 * @brief The C interface of sablib smoothing functions.
 * @author Izadori
 */

#ifndef __SABLIB_C_SMOOTHING_H__
#define __SABLIB_C_SMOOTHING_H__

#include "data.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * @brief Calculates the simple moving average of the input signal.
 *
 * @param y The data to be averaged.
 * @param n Half-width of the moving average window (calculated using `2 * n + 1` points).
 * @return The data after applying the moving average.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_MovingAverage(const SABLIB_DATA_PTR y, const unsigned int n);

/**
 * @brief Calculates the weighted moving average of the input signal.
 *
 * @param y The data to be averaged.
 * @param w Weights.
 * @return The data after applying the moving average.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_WeightedMovingAverage(const SABLIB_DATA_PTR y, const SABLIB_DATA_PTR w);

/**
 * @brief Generates a Gaussian kernel.
 *
 * @param n Half-width of the Gaussian window (total size is `2 * n + 1`).
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return A vector containing the Gaussian kernel coefficients.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_GaussianKernel(const unsigned int n, const double sigma);

/**
 * @brief Performs Gaussian smoothing on the input signal.
 *
 * This is a convenience function that generates a Gaussian kernel and then
 * applies it using `Sablib_WeightedMovingAverage`.
 *
 * @param y The data to be smoothed.
 * @param n Half-width of the Gaussian window (total size is `2 * n + 1`).
 * @param sigma The standard deviation of the Gaussian distribution.
 * @return The data after applying the Gaussian filter.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_GaussianFilter(const SABLIB_DATA_PTR y, const unsigned int n, const double sigma);

/**
 * @brief Performs moving median smoothing.
 *
 * @param y The input data vector (signal to be smoothed).
 * @param n Half-width of the moving median window (calculated using `2 * n + 1` points).
 * @return The data after applying the moving median.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_MovingMedian(const SABLIB_DATA_PTR y, const unsigned int n);

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
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_PSpline(
	const SABLIB_DATA_PTR y, const unsigned int knots_num,
	const unsigned int degree, const unsigned int s, const double lambda
);

/**
 * @brief Calculates the coefficients for a Savitzky-Golay filter.
 *
 * @param n The half-width of the filter window (total window size is `2 * n + 1`).
 * @param polyorder The order of the polynomial used for fitting.
 * @param derive The order of the derivative to compute (0 for smoothing only).
 * @param delta The spacing of the samples.
 * @return The calculated filter coefficients.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_SavitzkyGolayCoefficients(
	const unsigned int n, const unsigned int polyorder, const unsigned derive, const double delta
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
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_SavitzkyGolay(
	const SABLIB_DATA_PTR y,
	const unsigned int n, const unsigned int polyorder, const unsigned derive, const double delta
);

/**
 * @brief Performs Whittaker smoothing (with weights).
 *
 * @param y The input data to be smoothed.
 * @param w Weights for each data point.
 * @param lambda Smoothing parameter (larger values lead to more smoothing, but may flatten peaks).
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @return The smoothed data.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_WeightedWhittaker(
	const SABLIB_DATA_PTR y, const SABLIB_DATA_PTR w,
	const double lambda, const unsigned int s
);

/**
 * @brief Performs Whittaker smoothing (without weights).
 *
 * @param y The input data to be smoothed.
 * @param lambda Smoothing parameter (larger values lead to more smoothing, but may flatten peaks).
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @return The smoothed data.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_Whittaker(
	const SABLIB_DATA_PTR y, const double lambda, const unsigned int s
);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __SABLIB_C_SMOOTHING_H__
