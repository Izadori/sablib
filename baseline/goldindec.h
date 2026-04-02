/**
 * @file goldindec.h
 * @brief Baseline estimation using the Goldindec algorithm
 * @author Izadori
 * @details
 * References:
 * @li Liu, J.; Sun, J.; Huang, X.; Li, G.; Liu, B. "Goldindec: A Novel Algorithm for Raman Spectrum Baseline Correction" Appl. Spectrosc., 2015, 69(7), 834-842.
 */

#ifndef __SABLIB_GOLDINDEC_H__
#define __SABLIB_GOLDINDEC_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Performs baseline estimation using the Goldindec algorithm.
 *
 * @param y The input data vector (signal to be processed).
 * @param polyorder The order of the polynomial to be fitted.
 * @param peak_ratio Estimated ratio related to the peak content (default: 0.5).
 * @param alpha A weighting parameter for the iterative process (default: 0.99 * 0.5).
 * @param loop Maximum number of iterations for the main loop (default: 100).
 * @param eps Convergence threshold for the main loop (default: 1.0e-4).
 * @param loop_legend Maximum number of iterations for internal LEGEND algorithm (default: 50).
 * @param eps_legend Convergence threshold for internal LEGEND algorithm (default: 1e-3).
 * @param eps_s Convergence threshold for the internal parameter s (default: 1e-4).
 * @return A vector of the same size as y containing the estimated baseline.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double>
BaselineGoldindec(
	const std::vector<double> & y, const unsigned int polyorder, const double peak_ratio = 0.5,
	const double alpha = 0.99 * 0.5, const unsigned int loop = 100, const double eps = 1.0e-4,
	const unsigned int loop_legend = 50, const double eps_legend = 1e-3, const double eps_s = 1e-4
);

}; // namespace sablib

#endif // __SABLIB_GOLDINDEC_H__
