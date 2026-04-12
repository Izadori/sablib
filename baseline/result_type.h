/**
 * @file result_type.h
 * @brief Common definition of return value for baseline estimation functions
 * @author Izadori
 */

#ifndef __SABLIB_BASELINE_RESULT_TYPE_H__
#define __SABLIB_BASELINE_RESULT_TYPE_H__

#include <vector>

namespace sablib {

/**
 * @brief Return value structure for baseline estimation functions.
 */
struct BaselineResult
{
    std::vector<double> baseline;  /**< Baseline data. */
    std::vector<double> corrected; /**< Corrected(baseline-subtracted) data. */
};

}; // namespace sablib

#endif // __SABLIB_BASELINE_RESULT_TYPE_H__
