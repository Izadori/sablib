/**
 * @file aispls.cpp
 * @brief Baseline estimation using adaptive iterative smoothing parameter Penalized Least Squares(aisPLS)(implementation)
 * @author Izadori
 */

#include <cmath>
#include <functional>

#include "../misc/functions.h"
#include "../smoothing/whittaker.h"

#include "aispls.h"

namespace sablib {

//
// Implementation of BaselineAisPLS() function
//
const BaselineResult BaselineAisPLS(
	const std::vector<double> & y, const double lambda, const double r,
	const unsigned int s, const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineAisPLS(): the length of y is zero.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("BaselineAisPLS(): non-positive lambda value is given.");
	}

	if(r <= 0) {
		throw std::invalid_argument("BaselineAisPLS(): non-positive r value is given.");
	}

	if(s == 0 || s > 3) {
		throw std::invalid_argument("BaselineAisPLS(): s must be 1, 2 or 3.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineAisPLS(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineAisPLS(): non-positive eps value is given.");
	}

	size_t m = y.size();
	Eigen::VectorXd yy, w, l, wt, z, d, tmp1, tmp2;
	Eigen::SparseMatrix<double> I, D, L, DTD, D2;
	double d2y_norm, q, mean, sd;
	Eigen::Index ct;

	yy = Eigen::VectorXd::Map(y.data(), m);

	w.setOnes(m);
	l = Eigen::VectorXd::Constant(m, lambda);

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	DTD = D.transpose() * D;

	D2 = Diff(I, 2);
	d2y_norm = (D2 * yy).norm();

	for(unsigned int i = 0; i < loop; i++) {
		L = l.asDiagonal();
		z = Whittaker(yy, w, L * DTD);

		q = (D2 * z).norm() / d2y_norm;

		d = yy - z;

		ct = (d.array() < 0).count();
		mean = (d.array() < 0).select(d, 0).sum() / ct;
		sd = std::sqrt((d.array() < 0).select((d.array() - mean).matrix(), 0).squaredNorm() / (ct - 1));

		auto eliminated_d =	(d.array() < 0 && (mean - 3 * sd) < d.array() && d.array() < (mean + 3 * sd));
		ct = eliminated_d.count();
		mean = eliminated_d.select(d, 0).sum() / ct;
		sd = std::sqrt(eliminated_d.select((d.array() - mean).matrix(), 0).squaredNorm() / (ct -1));

		tmp1 = (-2 * (i + 1) * (d.array() + mean - 2 * sd) / sd).matrix().unaryExpr(std::ref(logistic));
		tmp2 = (-2 * (i + 1) * (d.array() + mean - 2 * sd) / q / sd).matrix().unaryExpr(std::ref(logistic));
		wt = (d.array() >= -mean + 2 * sd).select(tmp1, (d.array() < 0).select(1, tmp2));

		if ((w - wt).norm() / w.norm() < eps) {
			break;
		}

		l = (d.array() >= -mean + 2 * sd).select(l * r, l / r);
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
