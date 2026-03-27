/**
 * @file moving_median.h
 * @brief Smoothing using moving median
 * @author Izadori
 */

#ifndef __SABLIB_MOVING_MEDIAN_H__
#define __SABLIB_MOVING_MEDIAN_H__

#include <stdexcept>
#include <vector>

namespace sablib {

/**
 * @brief Performs moving median smoothing.
 *
 * @param y The input data vector (signal to be smoothed).
 * @param n Half-width of the moving median window (calculated using `2 * n + 1` points).
 * @return The data after applying the moving median.
 * @exception std::invalid_argument If the size of y or n is zero.
 */
const std::vector<double> MovingMedian(const std::vector<double> & y, const unsigned int n);

}; // namespace sablib

#endif // __SABLIB_MOVING_MEDIAN_H__
