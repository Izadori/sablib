/**
 * @file kajfoszkwiatek.cpp
 * @brief Baseline estimation using Kajfosz-Kwiatek method(implementation)
 * @author Izadori
 */

#include <limits>
#include <vector>

#include <Eigen/Dense>

#include "kajfoszkwiatek.h"

namespace sablib {

//
// Implementation of BaselineKajfoszKwiatek() function
//
const BaselineResult BaselineKajfoszKwiatek(const std::vector<double> & y, const double width)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineKajfoszKwiatek(): the length of y is zero.");
	}

	if(width <= 0) {
		throw std::invalid_argument("BaselineKajfoszKwiatek(): the width parameter must be positive.");
	}

	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd baseline(yy.size());

	const double inv = 1.0 / (width * width);

	std::vector<int> v(yy.size());
	std::vector<double> z(yy.size() + 1);

	int k = 0;

	v[0] = 0;
	z[0] = -std::numeric_limits<double>::infinity();
	z[1] =  std::numeric_limits<double>::infinity();

	// Lower envelope
	for (int q = 1; q < yy.size(); q++) {
		double s;

		while (true) {
			int p = v[k];

			double fq = yy(q);
			double fp = yy(p);
			double dp = p;
			double dq = q;

			s = ((fq - fp) / inv + dq * dq - dp * dp) / (2.0 * (dq - dp));

			if (s > z[k]) {
				break;
			}

			if (k == 0) {
				break;
			}

			k--;
		}

		if (s <= z[k]) {
			s = -std::numeric_limits<double>::infinity();
		}

		k++;

		v[k] = q;
		z[k] = s;
		z[k + 1] = std::numeric_limits<double>::infinity();
	}

	k = 0;
	for (int x = 0; x < yy.size(); x++) {
		while (z[k + 1] < x) {
			k++;
		}

		int p = v[k];
		double dx = double(x - p);

		baseline(x) = yy(p) + inv * dx * dx;
	}

	BaselineResult result;
	result.baseline.resize(baseline.size());
	result.corrected.resize(baseline.size());

	Eigen::VectorXd::Map(result.baseline.data(), baseline.size()) = baseline;
	Eigen::VectorXd::Map(result.corrected.data(), baseline.size()) = yy - baseline;

	return result;
}

}; // namespace sablib
