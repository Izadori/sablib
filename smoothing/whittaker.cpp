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
	if(y.size() == 0) {
		throw std::invalid_argument("Whittaker(): length of y is zero.");
	}

	if(y.size() != w.size()) {
		throw std::invalid_argument("Whittaker(): the sizes of y and w do not match.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("Whittaker(): non-positive lambda value is given.");
	}

	if(s > 3 || s == 0) {
		throw std::invalid_argument("Whittaker(): s must be 1, 2, or 3.");
	}

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
