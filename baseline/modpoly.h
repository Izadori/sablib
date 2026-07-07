/**
 * @file modpoly.h
 * @brief Baseline estimation using Modified Polynomial(ModPoly) method
 * @author Izadori
 * @details
 * References:
 * @li Lieber, C.A.; Mahadevan-Jansen, A. "Automated Method for Subtraction of Fluorescence from Biological Raman Spectra" Applied Spectroscopy 2003, 57(11), 1363-1367.
 */

#ifndef __SABLIB_MODPOLY_H__
#define __SABLIB_MODPOLY_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Estimates the baseline using the Modified Polynomial (ModPoly) method.
 *
 * @param y The input data points to be processed.
 * @param polyorder The order of the polynomial to fit.
 * @param loop The maximum number of iterations (default is 50).
 * @param eps The convergence tolerance (default is 1.0e-3).
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselineModPoly(
	const std::vector<double> & y, const unsigned int polyorder,
	const unsigned int loop = 50, const double eps = 1.0e-3
);

}; // namespace sablib

#endif // __SABLIB_MODPOLY_H__
