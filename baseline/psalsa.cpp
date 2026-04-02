/**
 * @file psalsa.cpp
 * @brief Baseline estimation using Peaked Signal’s Asymmetric Least Squares Algorithm(psalsa)(implementation)
 * @author Izadori
 */

#include "../smoothing/whittaker.h"

#include "psalsa.h"

namespace sablib {

//
// Implementation of BaselinePsalsa() function
//
const std::vector<double> BaselinePsalsa(
	std::vector<double> & y, const double lambda, const double p, const double k,
	const unsigned int s, const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselinePsalsa(): the length of y is zero.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("BaselinePsalsa(): non-positive lambda value is given.");
	}

	if(p <= 0) {
		throw std::invalid_argument("BaselinePsalsa(): non-positive p value is given.");
	}

	if(k <= 0) {
		throw std::invalid_argument("BaselinePsalsa(): non-positive k value is given.");
	}

	if(s == 0 || s > 3) {
		throw std::invalid_argument("BaselinePsalsa(): s must be 1, 2 or 3.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselinePsalsa(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselinePsalsa(): non-positive eps value is given.");
	}

	size_t m = y.size();
	Eigen::VectorXd yy, w, w_old, z, pv(m), npv;
	Eigen::SparseMatrix<double> I, D, lambdaDTD;

	yy = Eigen::VectorXd::Map(y.data(), m);

	w.setOnes(m);
	w_old.setZero(m);
	pv.fill(p);
	npv = (1 - pv.array()).matrix();

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	lambdaDTD = lambda * (D.transpose() * D);

	for(unsigned int i = 0; i < loop; i++) {
		z = Whittaker(yy, w, lambdaDTD);

		pv = (yy - z).unaryExpr([&](const double x) {
			return p * std::exp(-x / k);
		});
		w = (yy.array() > z.array()).select(pv, npv);

		if(((w.array() - w_old.array()).abs() < eps).all()) {
			break;
		}

		w_old = w;
	}

	std::vector<double> result(z.size());
	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;
}

}; // namespace sablib
