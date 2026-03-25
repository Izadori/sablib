/**
 * @file moving_average.cpp
 * @brief Smoothing using simple/weighted moving average(implementation)
 * @author Izadori
 */
#define _USE_MATH_DEFINES
#include <cmath>
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

//
// Implementation of WeightedMovingAverage() function, std::vector<double> version
//
const std::vector<double> WeightedMovingAverage(const std::vector<double> & y, const std::vector<double> & w)
{
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd ww = Eigen::VectorXd::Map(w.data(), w.size());

	Eigen::VectorXd result = WeightedMovingAverage(yy, ww);

	std::vector<double> result_y(result.size());
	Eigen::VectorXd::Map(result_y.data(), result_y.size()) = result;

	return result_y;
}

//
// Implementation of GaussianKernel() function
//
const std::vector<double> GaussianKernel(const unsigned int n, const double sigma)
{
	if(n == 0) {
		throw std::invalid_argument("GaussianKernel(): n is zero.");
	}

	if(sigma <= 0) {
		throw std::invalid_argument("GaussianKernel(): non-positive sigma value is given.");
	}

	int points = 2 * n + 1;
	std::vector<double> result(points);

	double coeff = 1.0 / (std::sqrt(2 * M_PI) * sigma);
	auto gaussian = [&](const double xx) {
		double t = xx / sigma;
		return coeff * std::exp(-t * t / 2);
	};

	double sum = 0;
	for(int i = -(int)n; i <= (int)n; i++) {
		result[i + n] = gaussian(i);
		sum += result[i + n];
	}

	for(int i = 0; i < points; i++) {
		result[i] /= sum;
	}

	return result;
}

}; // namespace sablib
