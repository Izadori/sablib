/**
 * @file whittaker.cpp
 * @brief Smoothing using Whittaker smoother(implementation)
 * @author Izadori
 */

#include "whittaker.h"

namespace sablib {

//
// Implementation of Whittaker() function, std::vector<double> and using vector w version
//
const std::vector<double> Whittaker(
	const std::vector<double> & y, const std::vector<double> & w,
	const double lambda, const unsigned int s
)
{
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd ww = Eigen::VectorXd::Map(w.data(), w.size());
	Eigen::VectorXd z;
	Eigen::SparseMatrix<double> I, D, lambdaDTD;

	size_t m = y.size();

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	lambdaDTD = lambda * (D.transpose() * D);

	z = Whittaker(yy, ww, lambdaDTD);

	std::vector<double> result(z.size());
	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;
}

//
// Implementation of Whittaker() function, std::vector<double> version
//
const std::vector<double> Whittaker(
	const std::vector<double> & y, const double lambda, const unsigned int s
)
{
	std::vector<double> w(y.size());

	for(size_t i = 0; i < w.size(); i++) {
		w[i] = 1;
	}

	return Whittaker(y, w, lambda, s);
}

}; // namespace sablib
