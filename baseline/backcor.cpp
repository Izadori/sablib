/**
 * @file backcor.cpp
 * @brief Baseline estimation using iterative polynomial fitting with a non-quadratic cost function(Backcor algorithm)(implementation)
 * @author Izadori
 */

#include <cmath>
#include <functional>

#include "../misc/polyfit.h"

#include "backcor.h"

namespace sablib {

//
// Implementation of BaselineBackcor() function.
//
const std::vector<double>
BaselineBackcor(
	const std::vector<double> & y, const unsigned int polyorder, const BackcorFunc func,
	const double s, const double alpha, const unsigned int loop, const double eps
)
{
	double max_index = (double)y.size() - 1;
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd xx(y.size());

	for(unsigned int i = 0; i < y.size(); i++) {
		xx(i) = i / max_index;
	}

	Eigen::MatrixXd V = Vandermonde(xx, polyorder);
	Eigen::LDLT<Eigen::MatrixXd> ldltV = (V.transpose() * V).ldlt();

	std::function<double(double)> derivative;

	switch(func) {
	case BackcorFunc::Huber:
		derivative = [&](const double xx) {
			if(xx <= -s) {
				return -2 * s;
			}
			else if(s <= xx) {
				return 2 * s;
			}
			else {
				return 2 * xx;
			}
		};
		break;
	case BackcorFunc::AHuber:
		derivative = [&](const double xx) {
			if(xx < s) {
				return 2 * xx;
			}
			else {
				return 2 * s;
			}
		};
		break;
	case BackcorFunc::TQuad:
		derivative = [&](const double xx) {
			if(std::fabs(xx) < s) {
				return 2 * xx;
			}
			else {
				return 0.0;
			}
		};
		break;
	case BackcorFunc::ATQuad:
		derivative = [&](const double xx) {
			if(xx < s) {
				return 2 * xx;
			}
			else {
				return 0.0;
			}
		};
		break;
	case BackcorFunc::Indec:
		derivative = [&](const double xx) {
			if(s <= xx) {
				return -s * s * s / (2 * xx * xx);
			}
			else if(xx <= s) {
				return s * s * s / (2 * xx * xx);
			}
			else {
				return 2 * xx;
			}
		};
		break;
	case BackcorFunc::AIndec:
		derivative = [&](const double xx) {
			if(xx < s) {
				return 2 * xx;
			}
			else {
				return -s * s * s / (2 * xx * xx);
			}
		};
		break;
	}

	Eigen::VectorXd coeff = ldltV.solve(V.transpose() * yy);
	Eigen::VectorXd z = PolyVal(coeff, V);

	for(unsigned int i = 0; i < loop; i++) {
		Eigen::VectorXd z_old = z;
		Eigen::VectorXd e = yy - z_old;
		Eigen::VectorXd d = -e + alpha * e.unaryExpr(derivative);
		Eigen::VectorXd coeff_new = ldltV.solve(V.transpose() * (yy + d));

		z = PolyVal(coeff_new, V);

		if((z - z_old).norm() / z_old.norm() < eps) {
			break;
		}
	}

	std::vector<double> result(z.size());
	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;
}

}; // namespace sablib
