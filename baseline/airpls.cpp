/**
 * @file airpls.cpp
 * @brief Baseline estimation using adaptive iteratively reweighted Penalized Least Squares(airPLS)(implementation)
 * @author Izadori
 */

#include "airpls.h"

namespace sablib {

//
// Implementation of BaselineAirPLS() function
//
const std::vector<double> BaselineAirPLS(
	std::vector<double> & y, const double lambda, const unsigned int s,
	const unsigned int loop, const double eps
)
{
	size_t m = y.size();
	Eigen::VectorXd yy, w, z, d;
	Eigen::SparseMatrix<double> I, D, lambdaDTD;
	double y_abs_sum, d_sum_abs;

	yy = Eigen::VectorXd::Map(y.data(), m);

	w.setOnes(m);
	y_abs_sum = yy.array().abs().matrix().sum();

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	lambdaDTD = lambda * (D.transpose() * D);

	for(unsigned int i = 0; i < loop; i++) {
		z = Whittaker(yy, w, lambdaDTD);

		d = (yy.array() >= z.array()).select(0, yy - z);
		d_sum_abs = std::fabs(d.sum());

		if (d_sum_abs < eps * y_abs_sum) {
			break;
		}

		w = (yy.array() >= z.array()).select(0, ((loop * d.array().abs()) / d_sum_abs).exp());
		w(0) = w(w.size() - 1) = std::exp((loop * d.maxCoeff() / d_sum_abs));
	}

	std::vector<double> result(z.size());

	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;

}

}; // namespace sablib
