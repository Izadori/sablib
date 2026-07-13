/**
 * @file functions.cpp
 * @brief Mathematical functions for baseline estimation(implementation)
 * @author Izadori
 */

#include "functions.h"

namespace sablib {

//
// Implementation of logistic() function
//
const double logistic(const double x)
{
	if(x < 0) {
		double t = std::exp(x);
		return t / (t + 1);
	}
	else {
		return 1 / (1 + std::exp(-x));
	}
}

}; // namespace sablib
