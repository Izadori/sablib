/**
 * @file imodpoly.cpp
 * @brief Baseline estimation using Improved Modified Polynomial(IModPoly) method(implementation)
 * @author Izadori
 */

#include <cmath>
#include <limits>

#include "../misc/polyfit.h"

#include "imodpoly.h"

namespace sablib {

//
// Implementation of BaselineIModPoly() function
//
const std::vector<double> BaselineIModPoly(
	const std::vector<double> & y, const unsigned int polyorder, const double k, const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineIModPoly(): the length of y is zero.");
	}

	if(polyorder == 0) {
		throw std::invalid_argument("BaselineIModPoly(): polyorder is zero.");
	}

	if(k <= 0) {
		throw std::invalid_argument("BaselineIModPoly(): non-positive k value is given.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineIModPoly(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineIModPoly(): non-positive eps value is given.");
	}

	double max_index = y.size() - 1;
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd xx(y.size());

	for(unsigned int i = 0; i < y.size(); i++) {
		xx(i) = i / max_index;
	}

	Eigen::MatrixXd V = Vandermonde(xx, polyorder);
	Eigen::LDLT<Eigen::MatrixXd> ldltV = (V.transpose() * V).ldlt();
	Eigen::VectorXd y_old = yy;
	double sd_old = std::numeric_limits<double>::infinity();

	for(unsigned int i = 0; i < loop; i++) {
		Eigen::VectorXd y_new = PolyVal(ldltV.solve(V.transpose() * y_old), V);

		Eigen::VectorXd r = y_new - y_old;
		double mean = r.mean();
		double sd = std::sqrt((r.array() - mean).square().sum() / r.size());

		if(std::fabs(sd - sd_old) < eps) {
			break;
		}

		y_old = (y_old.array() > (y_new.array() + k * sd)).select(y_new, y_old);
		sd_old = sd;
	}

	std::vector<double> result(y.size());
	Eigen::VectorXd::Map(result.data(), result.size()) = y_old;

	return result;
}

}; // namespace sablib
