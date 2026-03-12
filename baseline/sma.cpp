//
// sma.cpp
//
// Copyright (c) 2026 Izadori
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include "sma.h"

namespace sablib {

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
