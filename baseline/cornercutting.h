/**
 * @file cornercutting.h
 * @brief Baseline estimation using Corner-Cutting method
 * @author Izadori
 * @details
 * References:
 * @li Liu, Y.J.; Zhou X.G.; Yu Y.D. "A concise iterative method with Bezier technique for baseline construction" Analyist 2015, 140(23), 7984–7996.
 */

#ifndef __SABLIB_CORNERCUTTING_H__
#define __SABLIB_CORNERCUTTING_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs baseline estimation using Corner-Cutting method.
 *
 * @param y The input data for baseline estimation.
 * @param loop Maximum number of iterations.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument One or more parameters are wrong.
 */
SABLIB_EXPORT const BaselineResult BaselineCornerCutting(const std::vector<double> & y, const unsigned int loop = 100);


}; // namespace sablib

#endif // __SABLIB_CORNERCUTTING_H__
