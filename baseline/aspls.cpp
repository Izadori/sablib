/**
 * @file aspls.cpp
 * @brief Baseline estimation using adaptive smoothness Penalized Least Squares(asPLS)(implementation)
 * @author Izadori
 */

#include <cmath>
#include <functional>

#include "../misc/functions.h"
#include "../smoothing/whittaker.h"

#include "aspls.h"

namespace sablib {

//
// Implementation of BaselineAsPLS() function
//
const BaselineResult BaselineAsPLS(
	const std::vector<double> & y, const double lambda, const double k,
	const unsigned int s, const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineAsPLS(): the length of y is zero.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("BaselineAsPLS(): non-positive lambda value is given.");
	}

	if(k <= 0) {
		throw std::invalid_argument("BaselineAsPLS(): k must be a positive value.");
	}

	if(s == 0 || s > 3) {
		throw std::invalid_argument("BaselineAsPLS(): s must be 1, 2 or 3.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineAsPLS(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineAsPLS(): non-positive eps value is given.");
	}

	size_t m = y.size();
	Eigen::VectorXd yy, w, wt, alpha, z, d;
	Eigen::SparseMatrix<double> A, I, D, lambdaDTD;
	double mean, sd;
	Eigen::Index ct;

	yy = Eigen::VectorXd::Map(y.data(), m);

	w.setOnes(m);
	alpha.setOnes(m);

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	lambdaDTD = lambda * (D.transpose() * D);

	for(unsigned int i = 0; i < loop; i++) {
		A = alpha.asDiagonal();
		z = Whittaker(yy, w, A * lambdaDTD);

		d = yy - z;
		alpha = d.cwiseAbs();
		alpha /= alpha.maxCoeff();

		ct = (d.array() < 0).count();
		mean = (d.array() < 0).select(d, 0).sum() / ct;
		sd = std::sqrt((d.array() < 0).select((d.array() - mean).matrix(), 0).squaredNorm() / (ct - 1));

		wt = (-k * (d.array() - sd) / sd).matrix().unaryExpr(std::ref(logistic));

		if ((w - wt).norm() / w.norm() < eps) {
			break;
		}

		w = wt;
	}

	BaselineResult result;
	result.baseline.resize(z.size());
	result.corrected.resize(z.size());

	Eigen::VectorXd::Map(result.baseline.data(), z.size()) = z;
	Eigen::VectorXd::Map(result.corrected.data(), z.size()) = yy - z;

	return result;
}

}; // namespace sablib
