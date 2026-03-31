/**
 * @file backcor.h
 * @brief Baseline estimation using iterative polynomial fitting with a non-quadratic cost function(Backcor algorithm)
 * @author Izadori
 * @details
 * References:
 * @li Mazet, V.; Carteret, C.; Brie, D.; Idier, J.; Humbert, B. "Background removal from spectra by designing and minimising a non-quadratic cost function" Chemometrics and Intelligent Laboratory Systems, 2005, 76(2), 121-133.
 * @li Liu, J.; Sun, J.; Huang, X.; Li, G.; Liu, B. "Goldindec: A Novel Algorithm for Raman Spectrum Baseline Correction" Appl. Spectrosc., 2015, 69(7), 834-842.
 */

#ifndef __SABLIB_BACKCOR_H__
#define __SABLIB_BACKCOR_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Cost function types for the Backcor algorithm.
 *
 * These functions define how the fitting error is penalized during the iterative polynomial fitting process.
 * Asymmetric versions are typically used for signals with positive peaks to ensure the baseline stays below the signal.
 */
enum class BackcorFunc
{
	Huber,  /**< Huber function. */
	AHuber, /**< Asymmetric Huber function. */
	TQuad,  /**< Truncated quadratic function. */
	ATQuad, /**< Asymmetric truncated quadratic function. */
	Indec,  /**< Indec function. */
	AIndec  /**< Asymmetric Indec function. */
};

/**
 * @brief Performs baseline estimation using iterative polynomial fitting with a non-quadratic cost function (Backcor).
 *
 * This algorithm fits a polynomial of order polyorder to the signal y. In each iteration,
 * the signal is modified based on the chosen cost function func to progressively estimate the background.
 *
 * @param y The input data vector (signal to be processed).
 * @param polyorder The order of the polynomial to be fitted.
 * @param func The type of cost function to use (default: ATQuad).
 * @param s Threshold parameter for the cost function (default: 1.0).
 * @param alpha Control parameter for the iterative update (default: 0.99). Should be in range [0, 1].
 * @param loop The maximum number of iterations (default: 50).
 * @param eps Convergence threshold for the relative change in the estimated baseline (default: 1.0e-3).
 * @return A vector of the same size as y containing the estimated baseline.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double>
BaselineBackcor(
	const std::vector<double> & y, const unsigned int polyorder, const BackcorFunc func = BackcorFunc::ATQuad,
	const double s = 1, const double alpha = 0.99, const unsigned int loop = 50, const double eps = 1.0e-3
);

}; // namespace sablib

#endif // __SABLIB_BACKCOR_H__
