/**
 * @file pspline.cpp
 * @brief Smoothing with penalized Spline(P-Spline)(implementation)
 * @author Izadori
 */

#include "../misc/diff.h"
#include "../misc/bspline.h"

#include "pspline.h"

namespace sablib {

//
// Implementation of PSpline() function
//
const std::vector<double> PSpline(
	const std::vector<double> & y, const unsigned int knots_num,
	const unsigned int degree, const unsigned int s, const double lambda
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("PSpline(): the length of y is zero.");
	}

	if(knots_num == 0) {
		throw std::invalid_argument("PSpline(): knots_num is zero.");
	}

	if(degree == 0) {
		throw std::invalid_argument("PSpline(): degree is zero.");
	}

	if(y.size() < knots_num) {
		throw std::invalid_argument("PSpline(): knots_num is larger than length of y.");
	}

	if(s == 0 || s > 3) {
		throw std::invalid_argument("PSpline(): s must be 1, 2 or 3.");
	}

	if(lambda <= 0) {
		throw std::invalid_argument("PSpline(): non-positive lambda value is given.");
	}

	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd xx = Eigen::VectorXd::LinSpaced(yy.size(), 0, 1);

	Eigen::VectorXd interior_knots = Eigen::VectorXd::LinSpaced(knots_num, 0, 1);
	Eigen::VectorXd knots(knots_num + degree * 2);

	for(unsigned int i = 0; i < degree; i++) {
		knots(i) = interior_knots(0);
		knots(knots.size() - 1 - i) = interior_knots(knots_num - 1);
	}

	knots.block(degree, 0, knots_num, 1) = interior_knots;

	BSpline bs(degree, knots);
	Eigen::SparseMatrix<double> B = bs.DesignMatrix(xx);
	Eigen::SparseMatrix<double> I, D, lambdaDTD, BTB;

	I.resize(B.cols(), B.cols());
	I.setIdentity();
	D = Diff(I, s);
	lambdaDTD = lambda * (D.transpose() * D);
	BTB = B.transpose() * B;

	Eigen::SimplicialCholesky< Eigen::SparseMatrix<double> > solver;

	solver.compute(BTB + lambdaDTD);

	if(solver.info() != Eigen::Success) {
		throw std::runtime_error("PSpline(): solver calculation fails.");
	}

	Eigen::VectorXd coefficients = solver.solve(B.transpose() * yy);
	std::vector<double> result(yy.size());

	for(int i = 0; i < xx.size(); i++) {
		result[i] = bs.Interpolate(xx(i), coefficients);
	}

	return result;
}

}; // namespace sablib
