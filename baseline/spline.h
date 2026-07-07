/**
 * @file spline.h
 * @brief Baseline estimation from cubic spline
 * @author Izadori
 */

#ifndef __SABLIB_SPLINE_H__
#define __SABLIB_SPLINE_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using cubic spline interpolation.
 *
 * @param y The input data for baseline estimation.
 * @param indices The indices of the points to be used as knots for the cubic spline.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselineSpline(const std::vector<double> & y, const std::vector<unsigned int> & indices);

}; // namespace sablib

#endif // __SABLIB_SPLINE_H__
