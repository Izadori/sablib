/**
 * @file arpls.cpp
 * @brief Baseline estimation using asymmetrically reweighted Penalized Least Squares(arPLS)(implementation)
 * @author Izadori
 */

#include <cmath>
#include <limits>

#include "arpls.h"

namespace sablib {

//
// Implementation of BaselineArPLS() function
//
const std::vector<double> BaselineArPLS(
	std::vector<double> & y, const double lambda, const unsigned int s,
	const unsigned int loop, const double eps
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineArPLS(): the length of y is zero.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("BaselineArPLS(): non-positive lambda value is given.");
	}

	if(s == 0 || s > 3) {
		throw std::invalid_argument("BaselineArPLS(): s must be 1, 2 or 3.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineArPLS(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineArPLS(): non-positive eps value is given.");
	}

	size_t m = y.size();
	Eigen::VectorXd yy, w, wt, z, d, tmp;
	Eigen::SparseMatrix<double> I, D, lambdaDTD;
	double mean, sd;
	Eigen::Index ct;

	const double overflow = std::log(std::numeric_limits<double>::max()) / 2;

	yy = Eigen::VectorXd::Map(y.data(), m);

	w.setOnes(m);

	I.resize(m, m);
	I.setIdentity();
	D = Diff(I, s);
	lambdaDTD = lambda * (D.transpose() * D);

	for(unsigned int i = 0; i < loop; i++) {
		z = Whittaker(yy, w, lambdaDTD);

		d = yy - z;

		ct = (d.array() < 0).count();
		mean = (d.array() < 0).select(d, 0).sum() / ct;
		sd = std::sqrt((d.array() < 0).select((d.array() - mean).matrix(), 0).squaredNorm() / (ct - 1));

		tmp = ((d.array() + mean - 2 * sd) / sd).matrix();
		tmp = (tmp.array() >= overflow).select(0, (1.0 / (1 + (2 * tmp).array().exp())).matrix());

		wt = (d.array() >= 0).select(tmp, 1);

		if ((w - wt).norm() / w.norm() < eps) {
			break;
		}

		w = wt;
	}

	std::vector<double> result(z.size());

	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;
}

}; // namespace sablib
