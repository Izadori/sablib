/**
 * @file sma.cpp
 * @brief Baseline estimation using simple moving average(implementation)
 * @author Izadori
 */

#include "sma.h"

namespace sablib {

//
// Implementation of BaselineSMA() function
//
std::vector<double> BaselineSMA(std::vector<double> & y, const unsigned int n, const unsigned int loop)
{
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd result = yy;

	for(unsigned int i = 0; i < loop; i++){
		Eigen::VectorXd result_old = result;
		result = MovingAverage(result, n);
		result = (result.array() > result_old.array()).select(result_old, result);
	}

	std::vector<double> result_y(result.size());
	Eigen::VectorXd::Map(result_y.data(), result_y.size()) = result;

	return result_y;
}

}; // namespace sablib
