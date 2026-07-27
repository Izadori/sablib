/**
 * @file rubberband.h
 * @brief Baseline estimation using Rubber Band method
 * @author Izadori
 */

#ifndef __SABLIB_RUBBERBAND_H__
#define __SABLIB_RUBBERBAND_H__

#include <stdexcept>

#include "result_type.h"
#include "sablib_export.h"

namespace sablib {

/**
 * @brief Performs background estimation using the Rubber Band method.
 *
 * @param y The input data for baseline estimation.
 * @return The estimated baseline and baseline-corrected data.
 * @exception std::invalid_argument Thrown if the input vector is empty.
 */
SABLIB_EXPORT const BaselineResult BaselineRubberBand(const std::vector<double> & y);

}; // namespace sablib

#endif // __SABLIB_RUBBERBAND_H__
