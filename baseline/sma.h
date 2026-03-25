/**
 * @file sma.h
 * @brief Baseline estimation using simple moving average
 * @author Izadori
 */

#ifndef __SABLIB_SMA_H__
#define __SABLIB_SMA_H__

#include <vector>

namespace sablib {

/**
 * @brief Performs background estimation using a simple moving average.
 *
 * @param y The input data for baseline estimation.
 * @param n Half-width of the moving average window (calculated using `2 * n + 1` points).
 * @param loop Number of iterations (default is 50).
 * @return The estimated baseline.
 */
const std::vector<double> BaselineSMA(std::vector<double> & y, const unsigned int n, const unsigned int loop = 50);

}; // namespace sablib

#endif // __SABLIB_SMA_H__
