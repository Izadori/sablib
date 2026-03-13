/**
 * @file moving_average.cpp
 * @brief Smoothing using simple moving average(implementation)
 * @author Izadori
 */

#include "moving_average.h"

namespace sablib {

//
// Implementation of MovingAverage() function, std::vector<double> version
//
const std::vector<double> MovingAverage(const std::vector<double> & y, const unsigned int n)
{
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());

	Eigen::VectorXd result = MovingAverage(yy, n);

	std::vector<double> result_y(result.size());
	Eigen::VectorXd::Map(result_y.data(), result_y.size()) = result;

	return result_y;
}

}; // namespace sablib
