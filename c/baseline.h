/**
 * @file baseline.h
 * @brief The C interface of sablib baseline estimation functions.
 * @author Izadori
 */

#ifndef __SABLIB_C_BASELINE_H__
#define __SABLIB_C_BASELINE_H__

#ifndef __cplusplus
    #include <stdbool.h>
#endif

#include "data.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * @brief Cost function types for the Backcor algorithm.
 */
enum Sablib_BackcorFunc
{
    Huber = 0,  /**< Huber function. */
    AHuber,     /**< Asymmetric Huber function. */
    TQuad,      /**< Truncated quadratic function. */
    ATQuad,     /**< Asymmetric truncated quadratic function. */
    Indec,      /**< Indec function. */
    AIndec      /**< Asymmetric Indec function. */
};

/**
 * @brief Preprocessing types for the SNIP algorithm.
 */
enum Sablib_SnipPreprocess
{
    None = 0,   /**< No preprocessing (linear scale). */
    LL,         /**< Log-log transformation. */
    LLS         /**< Log-log-sqrt transformation. */
};

/**
 * @brief Penalty types for the BEADS algorithm.
 */
enum Sablib_BeadsPenalty
{
    L1_v1 = 0,  /**< Penalty type 1. */
    L1_v2       /**< Penalty type 2. */
};

/**
 * @brief Performs baseline estimation with a linear line between two points.
 *
 * @param y The input data for baseline estimation.
 * @param index1 The index of the first point.
 * @param index2 The index of the second point.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineLinear(
    const SABLIB_DATA_PTR y, const unsigned int index1, const unsigned int index2
);

/**
 * @brief Performs baseline estimation by fitting a polynomial to specified points.
 *
 * @param y The input data for baseline estimation.
 * @param polyorder The order of the polynomial to fit.
 * @param indices_ptr The indices of the points used for polynomial fitting.
 * @param ptr_size The number of indices.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselinePolynomial(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const unsigned int * indices_ptr, const size_t ptr_size
);

/**
 * @brief Performs baseline estimation using cubic spline interpolation.
 *
 * @param y The input data for baseline estimation.
 * @param indices_ptr The indices of the points to be used as knots for the cubic spline.
 * @param ptr_size The number of indices.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineSpline(
    const SABLIB_DATA_PTR y, const unsigned int * indices_ptr, const size_t ptr_size
);

/**
 * @brief Performs background estimation using a simple moving average.
 *
 * @param y The input data for baseline estimation.
 * @param n Half-width of the moving average window (calculated using `2 * n + 1` points).
 * @param loop Number of iterations.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineSMA(
    const SABLIB_DATA_PTR y, const unsigned int n, const unsigned int loop
);

/**
 * @brief Performs baseline estimation using the Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP) algorithm.
 *
 * @param y The input data vector (signal to be processed).
 * @param m The maximum half-window size (maximum clipping distance).
 * @param decreasing If true, iterates from m down to 1 (recommended). If false, iterates from 1 up to m.
 * @param preprocess The preprocessing transformation to apply before clipping.
 * @param loop The number of times to repeat the entire SNIP process.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineSnip(
	const SABLIB_DATA_PTR y, const unsigned int m, const bool decreasing,
	const enum Sablib_SnipPreprocess preprocess, const unsigned int loop
);

/**
 * @brief Estimates the baseline using the Modified Polynomial (ModPoly) method.
 *
 * @param y The input data points to be processed.
 * @param polyorder The order of the polynomial to fit.
 * @param loop The maximum number of iterations.
 * @param eps The convergence tolerance.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineModPoly(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const unsigned int loop, const double eps
);

/**
 * @brief Estimates the baseline using the Improved Modified Polynomial (IModPoly) method.
 *
 * @param y The input data points to be processed.
 * @param polyorder The order of the polynomial to fit.
 * @param k Scaling factor for the standard deviation threshold.
 * @param loop The maximum number of iterations.
 * @param eps The convergence tolerance based on the standard deviation change.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineIModPoly(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const double k,
	const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using iterative polynomial fitting with a non-quadratic cost function (Backcor).
 *
 * @param y The input data vector (signal to be processed).
 * @param polyorder The order of the polynomial to be fitted.
 * @param func The type of cost function to use.
 * @param s Threshold parameter for the cost function.
 * @param alpha Control parameter for the iterative update. Should be in range [0, 1].
 * @param loop The maximum number of iterations.
 * @param eps Convergence threshold for the relative change in the estimated baseline.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineBackcor(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const enum Sablib_BackcorFunc func,
	const double s, const double alpha, const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using the Goldindec algorithm.
 *
 * @param y The input data vector (signal to be processed).
 * @param polyorder The order of the polynomial to be fitted.
 * @param peak_ratio Estimated ratio related to the peak content.
 * @param alpha A weighting parameter for the iterative process.
 * @param loop Maximum number of iterations for the main loop.
 * @param eps Convergence threshold for the main loop.
 * @param loop_legend Maximum number of iterations for internal LEGEND algorithm.
 * @param eps_legend Convergence threshold for internal LEGEND algorithm.
 * @param eps_s Convergence threshold for the internal parameter s.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineGoldindec(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const double peak_ratio,
	const double alpha, const unsigned int loop, const double eps,
	const unsigned int loop_legend, const double eps_legend, const double eps_s
);

/**
 * @brief Performs baseline estimation using Asymmetric Least Squares Smoothing (AsLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param p Weight (asymmetry parameter, typically 0.001 to 0.1).
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAsLS(
	const SABLIB_DATA_PTR y, const double lambda, const double p, const unsigned int s,
	const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using Improved Asymmetric Least Squares Smoothing (IAsLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param lambda1 Smoothing parameter for the first derivative.
 * @param p Weight (asymmetry parameter, typically 0.001 to 0.1).
 * @param s The order of the difference (usually s = 2 or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineIAsLS(
	const SABLIB_DATA_PTR y, const double lambda, const double lambda1, const double p,
	const unsigned int s, const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using adaptive iteratively reweighted Penalized Least Squares(airPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAirPLS(
	const SABLIB_DATA_PTR y, const double lambda, const unsigned int s,
	const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using adaptive iterative smoothing parameter Penalized Least Squares(aisPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param r Adaptation rate of smoothing parameter.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAisPLS(
	const SABLIB_DATA_PTR y, const double lambda, const double r,
	const unsigned int s, const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using asymmetrically reweighted Penalized Least Squares(arPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineArPLS(
	const SABLIB_DATA_PTR y, const double lambda, const unsigned int s,
	const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using adaptive smoothness Penalized Least Squares(asPLS).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param k Asymmetric coefficient.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAsPLS(
	const SABLIB_DATA_PTR y, const double lambda, const double k,
	const unsigned int s, const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation using Peaked Signal’s Asymmetric Least Squares Algorithm (psalsa).
 *
 * @param y The input data for baseline estimation.
 * @param lambda Smoothing parameter.
 * @param p Weight (asymmetry parameter, typically 0.001 to 0.1).
 * @param k Exponential decay of the weights.
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselinePsalsa(
	const SABLIB_DATA_PTR y, const double lambda, const double p, const double k,
	const unsigned int s, const unsigned int loop, const double eps
);

/**
 * @brief Performs baseline estimation and denoising using Sparsity (BEADS).
 *
 * @param y The input data.
 * @param s Order of the derivative for baseline sparsity (typically 1 or 2).
 * @param frequency Sampling frequency of the signal.
 * @param r High-pass filter parameter (cut-off frequency relative to sampling frequency).
 * @param lambda0 Sparsity parameter for the baseline.
 * @param lambda1 Sparsity parameter for the first-order derivative of the signal.
 * @param lambda2 Sparsity parameter for the second-order derivative of the signal.
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @param penalty Penalty type.
 * @return The estimated baseline and denoised baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineBeads(
	const SABLIB_DATA_PTR y, const unsigned int s, const double frequency, const double r,
	const double lambda0, const double lambda1, const double lambda2, const unsigned int loop,
	const double eps, const enum Sablib_BeadsPenalty penalty
);

/**
 * @brief Expands the signal boundaries by padding with a tapered sequence.
 *
 * @param y The input data vector.
 * @param n The number of elements to add at each end.
 * @return A expanded vector.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_BeadsExpandBoundaries(const SABLIB_DATA_PTR y, const unsigned int n);

/**
 * @brief Trims the expanded signal boundaries.
 *
 * @param y The input data vector (expanded signal).
 * @param n The number of elements to trim from each end.
 * @return A trimmed vector.
 * @note The returned pointer must be freed with FreeSablibData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR Sablib_BeadsTrimBoundaries(const SABLIB_DATA_PTR y, const unsigned int n);

/**
 * @brief Performs background estimation using the Rubber Band method.
 *
 * @param y The input data for baseline estimation.
 * @return The estimated baseline and baseline-corrected data.
 * @note The returned pointer must be freed with FreeSablibBaselineData() to avoid memory leaks.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR Sablib_BaselineRubberBand(const SABLIB_DATA_PTR y);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __SABLIB_C_BASELINE_H__
