/**
 * @file kajfoszkwiatek.h
 * @brief Baseline estimation using Kajfosz-Kwiatek method
 * @author Izadori
 * References:
 * @li Kajfosz, J.; Kwiatek, W. "Nonpolynomial approximation of background in X-ray spectra" Nucl. Instrum. Methods 1987, B22, 78-81.
 * @li Felzenszwalb, P. F.; Huttenlocher, D. P. "Distance Transforms of Sampled Functions" Theory of Computing 2012, 19(8), 415-428.
 */

#ifndef __SABLIB_KAJFOSZ_KWIATEK_H__
#define __SABLIB_KAJFOSZ_KWIATEK_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs background estimation using the Kajfosz-Kwiatek method. (parabolic approximation)
 *
 * Using the Felzenszwalb-Huttenlocher's lower envelope algorithm.
 *
 * @param y The input data for baseline estimation.
 * @param width The width of parabola.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument Thrown if the input vector is empty.
 */
SABLIB_EXPORT const BaselineResult BaselineKajfoszKwiatek(const std::vector<double> & y, const double width);

}; // namespace sablib

#endif // __SABLIB_KAJFOSZ_KWIATEK_H__
