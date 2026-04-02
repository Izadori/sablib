/**
 * @file goldindec.cpp
 * @brief Baseline estimation using Goldindec(implementation)
 * @author Izadori
 */

#include <cmath>
#include <functional>
#include <iostream>

#include "../misc/polyfit.h"

#include "goldindec.h"

namespace sablib {

const std::vector<double>
BaselineGoldindec(
	const std::vector<double> & y, const unsigned int polyorder,
	double peak_ratio, const double alpha, const unsigned int loop, const double eps,
	const unsigned int loop_legend, const double eps_legend, const double eps_s
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineGoldindec(): the length of y is zero.");
	}

	if(polyorder == 0) {
		throw std::invalid_argument("BaselineGoldindec(): polyorder is zero.");
	}

	if(peak_ratio <= 0 || 1 < peak_ratio) {
		throw std::invalid_argument("BaselineGoldindec(): peak_ratio must be between 0 and 1.");
	}

	if(alpha <= 0 || 1 < alpha) {
		throw std::invalid_argument("BaselineGoldindec(): alpha must be between 0 and 1.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineGoldindec(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineGoldindec(): non-positive eps value is given.");
	}

	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(y.size(), 0, 1);
	Eigen::MatrixXd V = Vandermonde(x, polyorder);
	Eigen::LDLT<Eigen::MatrixXd> ldltV = (V.transpose() * V).ldlt();

	double a = 0, b = 1;
	double rud = 0.7679
		+ 11.2358 * peak_ratio
		- 39.7064 * peak_ratio * peak_ratio
		+ 92.3583 * peak_ratio * peak_ratio * peak_ratio; // See Liu's paper.
	double s = a + 0.618 * (b - a), s_old = s;

	auto derivative = [&](const double xx) {
		if(xx < s) {
			return 2 * xx;
		}
		else {
			return -s * s * s / (2 * xx * xx);
		}
	};

	Eigen::VectorXd coeff = ldltV.solve(V.transpose() * yy);
	Eigen::VectorXd z = PolyVal(coeff, V);

	for(unsigned int i = 0; i < loop; i++) {
		// LEGEND algorithm
		for(unsigned int j = 0; j < loop_legend; j++) {
			Eigen::VectorXd z_old = z;

			Eigen::VectorXd e = yy - z_old;
			Eigen::VectorXd d = -e + alpha * e.unaryExpr(derivative);
			coeff = ldltV.solve(V.transpose() * (yy + d));
			z = PolyVal(coeff, V);

			if((z - z_old).norm() / z_old.norm() < eps_legend) {
				break;
			}
		}

		double up_down_ratio = (double)(yy.array() >= z.array()).count() / (double)(yy.array() < z.array()).count();
		double diff = up_down_ratio - rud;

		if(diff > eps) {
			a = s;
		}
		else if(diff < -eps) {
			b = s;
		}
		else {
			break;
		}

		s = a + 0.618 * (b - a);
		// std::cerr << "up_down_ratio = " << up_down_ratio << ", s = " << s << std::endl;

		// Based on experimentation, it seems better to set a threshold for s as well.
		if(std::fabs(s - s_old) < eps_s) {
			break;
		}

		s_old = s;
	}

	std::vector<double> result(z.size());
	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;
}

}; // namespace sablib
