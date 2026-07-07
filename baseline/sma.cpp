/**
 * @file sma.cpp
 * @brief Baseline estimation using simple moving average(implementation)
 * @author Izadori
 */

#include "../smoothing/moving_average.h"

#include "sma.h"

namespace sablib {

//
// Implementation of BaselineSMA() function
//
const BaselineResult BaselineSMA(const std::vector<double> & y, const unsigned int n, const unsigned int loop)
{
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd result = yy;

	for(unsigned int i = 0; i < loop; i++){
		Eigen::VectorXd result_old = result;
		result = MovingAverage(result, n);
		result = (result.array() > result_old.array()).select(result_old, result);
	}

	BaselineResult result_y;
	result_y.baseline.resize(yy.size());
	result_y.corrected.resize(yy.size());

	Eigen::VectorXd::Map(result_y.baseline.data(), yy.size()) = result;
	Eigen::VectorXd::Map(result_y.corrected.data(), yy.size()) = yy - result;

	return result_y;
}

}; // namespace sablib
