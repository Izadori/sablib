/**
 * @file modpoly.cpp
 * @brief Baseline estimation using Modified Polynomial(ModPoly) method(implementation)
 * @author Izadori
 */

#include "../misc/polyfit.h"

#include "modpoly.h"

namespace sablib {

//
// Implementation of BaselineModPoly() function
//
const std::vector<double> BaselineModPoly(
	const std::vector<double> & y, const unsigned int polyorder, const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineModPoly(): the length of y is zero.");
	}

	if(polyorder == 0) {
		throw std::invalid_argument("BaselineModPoly(): polyorder is zero.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineModPoly(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineModPoly(): non-positive eps value is given.");
	}

	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(y.size(), 0, 1);
	Eigen::MatrixXd V = Vandermonde(x, polyorder);
	Eigen::LDLT<Eigen::MatrixXd> ldltV = (V.transpose() * V).ldlt();
	Eigen::VectorXd y_old = yy;

	for(unsigned int i = 0; i < loop; i++) {
		Eigen::VectorXd y_new = PolyVal(ldltV.solve(V.transpose() * y_old), V);

		if((y_new - y_old).norm() / y_old.norm() < eps) {
			break;
		}

		y_old = (y_new.array() > y_old.array()).select(y_old, y_new);
	}

	std::vector<double> result(y.size());
	Eigen::VectorXd::Map(result.data(), result.size()) = y_old;

	return result;
}

}; // namespace sablib
