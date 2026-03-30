/**
 * @file snip.cpp
 * @brief Baseline estimation using Statistics-sensitive Non-linear Iterative Peak-clipping(SNIP)(implementation)
 * @author Izadori
 */

#include <algorithm>
#include <cmath>

#include <Eigen/Eigen>

#include "snip.h"

namespace sablib {

const std::vector<double>
BaselineSnip(
	const std::vector<double> & y, const unsigned int m, const bool decreasing,
	const SnipPreprocess preprocess, const unsigned int loop
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineSnip(): the length of y is zero.");
	}

	if(y.size() < 2 * m) {
		throw std::invalid_argument("BaselineSnip(): the length of y is shorter than 2 * m.");
	}

	if(m == 0) {
		throw std::invalid_argument("BaselineSnip(): m is zero.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineSnip(): loop is zero.");
	}

	Eigen::VectorXd yy;
	unsigned int n = (unsigned int)y.size();

	switch(preprocess) {
	case SnipPreprocess::LL:
		yy = Eigen::VectorXd::Map(y.data(), n).unaryExpr(
			[](const double x) {
				return x > 0 ? std::log(std::log(x + 1) + 1) : 0;
			}
		);
		break;
	case SnipPreprocess::LLS:
		yy = Eigen::VectorXd::Map(y.data(), n).unaryExpr(
			[](const double x) {
				return x > 0 ? std::log(std::log(std::sqrt(x + 1) + 1) + 1) : 0;
			}
		);
		break;
	default:
		yy = Eigen::VectorXd::Map(y.data(), n);
	}

	std::vector<unsigned int> index;

	for(unsigned int i = 0; i < m; i++){
		index.emplace_back(decreasing ? (m - i) : (i + 1));
	}

	Eigen::VectorXd v = yy;

	for(unsigned int i = 0; i < loop; i++) {
		for(auto j : index) {
			Eigen::VectorXd v_prev = v;

			for(unsigned int k = j; k < n - j; k++) {
				double ave = (v_prev(k + j) + v_prev(k - j)) / 2;
				v(k) = std::min(v_prev(k), ave);
			}
		}
	}

	std::vector<double> result(v.size());

	switch(preprocess) {
	case SnipPreprocess::LL:
		Eigen::VectorXd::Map(result.data(), result.size()) = v.unaryExpr(
			[](const double x) {
				return std::exp(std::exp(x) - 1) - 1;
			}
		);
		break;
	case SnipPreprocess::LLS:
		Eigen::VectorXd::Map(result.data(), result.size()) = v.unaryExpr(
			[](const double x) {
				double t = std::exp(std::exp(x) - 1) - 1;
				return t * t  - 1;
			}
		);
		break;
	default:
		Eigen::VectorXd::Map(result.data(), result.size()) = v;
	}

	return result;
}

}; // namespace sablib
