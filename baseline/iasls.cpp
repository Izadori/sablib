/**
 * @file iasls.cpp
 * @brief Baseline estimation using Improved Asymmetric Least Squares Smoothing(IAsLS)(implementation)
 * @author Izadori
 */

#include "../misc/polyfit.h"
#include "../smoothing/whittaker.h"

#include "iasls.h"

namespace sablib {

//
// Implementation of BaselineIAsLS() function
//
const BaselineResult BaselineIAsLS(
	const std::vector<double> & y, const double lambda, const double lambda1, const double p,
	const unsigned int s, const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineIAsLS(): the length of y is zero.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("BaselineIAsLS(): non-positive lambda value is given.");
	}

	if(lambda1 <= 0) {
		throw std::invalid_argument("BaselineIAsLS(): non-positive lambda1 value is given.");
	}

	if(p <= 0) {
		throw std::invalid_argument("BaselineIAsLS(): non-positive p value is given.");
	}

	if(s <= 1 || s > 3) {
		throw std::invalid_argument("BaselineIAsLS(): s must be 2 or 3.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineIAsLS(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineIAsLS(): non-positive eps value is given.");
	}

	size_t m = y.size();
	Eigen::VectorXd yy, y_old, w, w1, z, z0, pv(m), npv;
	Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(y.size(), 0, 1);
	Eigen::SparseMatrix<double> I, D, D1, lambda1D1TD1, lambdaDTD, W, WTW;
	Eigen::MatrixXd V = Vandermonde(x, 2);
	Eigen::LDLT<Eigen::MatrixXd> ldltV = (V.transpose() * V).ldlt();
	Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;

	yy = Eigen::VectorXd::Map(y.data(), m);

	w.setOnes(m);
	z0.setZero(m);
	pv.fill(p);
	npv = (1 - pv.array()).matrix();

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	D1 = Diff(I, 1);
	lambda1D1TD1 = lambda1 * (D1.transpose() * D1);
	lambdaDTD = lambda * (D.transpose() * D);

	for(unsigned int i = 0; i < loop; i++) {
		z = PolyVal(ldltV.solve(V.transpose() * yy), V);
		w = (yy.array() > z.array()).select(pv, npv);

		y_old = yy;
		while(true) {
			W = w.asDiagonal();
			WTW = W.transpose() * W;

			solver.compute(WTW + lambdaDTD + lambda1D1TD1);
			z = solver.solve((WTW + lambda1D1TD1) * yy);

			if(solver.info() != Eigen::Success) {
				throw std::runtime_error("BaselineIAsLS(): solver calculation fails.");
			}

			w1 = (yy.array() > z.array()).select(pv, npv);

			if((w - w1).norm() / w.norm() < eps) {
				yy -= z;
				break;
			}

			w = w1;
		}

		z0 += z;

		if(((yy.array() - y_old.array()).abs() < eps).all()) {
			break;
		}
	}

	BaselineResult result;
	result.baseline.resize(z0.size());
	result.corrected.resize(z0.size());

	Eigen::VectorXd::Map(result.baseline.data(), z0.size()) = z0;
	Eigen::VectorXd::Map(result.corrected.data(), z0.size()) = yy - z0;

	return result;
}

}; // namespace sablib
