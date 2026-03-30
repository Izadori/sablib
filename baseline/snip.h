/**
 * @file snip.h
 * @brief Baseline estimation using Statistics-sensitive Non-linear Iterative Peak-clipping(SNIP)
 * @author Izadori
 * @details
 * References:
 * @li Ryan, C. G.; Clayton, E.; Griffin, W. L.; Sie, S. H.; Cousens, D. R. "SNIP, A Statistics-sensitive Background Treatment for the Quantitative Analysis of PIXE Spectra in Geoscience Applications" Nuclear Instruments and Methods in Physics Research, 1988, B34(3), 396-402.
 * @li Caccia, M.; Ebolese, A.; Matteo, M.; Santro, R.; Locatteli, M.; Pieracci, M.; Tintori, C. "Background removal procedure based on the SNIP algorithm for γ−ray spectroscopy with the CAEN Educational Kit" Educational Note ED3163; CAEN S.p.A.: Viareggio, Italy, 2021.
 */

#include <stdexcept>
#include <vector>

#ifndef __SABLIB_SNIP_H__
#define __SABLIB_SNIP_H__

namespace sablib {

/**
 * @brief Preprocessing types for the SNIP algorithm.
 *
 * SNIP often works better when the data is transformed to stabilize the variance
 * or handle high dynamic range peaks.
 */
enum class SnipPreprocess
{
	None, /**< No preprocessing (linear scale). */
	LL,   /**< Log-log transformation: `v = log(log(y + 1) + 1)`. Useful for large dynamic ranges. */
	LLS   /**< Log-log-sqrt transformation: `v = log(log(sqrt(y + 1) + 1) + 1)`. More aggressive variance stabilization. */
};

/**
 * @brief Performs baseline estimation using the Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP) algorithm.
 *
 * @param y The input data vector (signal to be processed).
 * @param m The maximum half-window size (maximum clipping distance).
 * @param decreasing If true, iterates from m down to 1 (recommended). If false, iterates from 1 up to m.
 * @param preprocess The preprocessing transformation to apply before clipping (None, LL, or LLS).
 * @param loop The number of times to repeat the entire SNIP process (usually 1 is sufficient).
 * @return A vector of the estimated baseline.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
const std::vector<double>
BaselineSnip(
	const std::vector<double> & y, const unsigned int m, const bool decreasing = true,
	const SnipPreprocess preprocess = SnipPreprocess::None, const unsigned int loop = 1
);

}; // namespace sablib

#endif // __SABLIB_SNIP_H__
