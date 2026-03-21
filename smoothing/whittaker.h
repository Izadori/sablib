/**
 * @file whittaker.h
 * @brief Smoothing using Whittaker smoother
 * @author Izadori
 * @details
 * References:
 * @li Whittaker, E. T. On A New Method of Graduation. Proc. Edinburgh Math. Soc., 1922, 41, 63-75.
 * @li Eilers, P. A Perfect Smoother. Analytical Chemistry, 2003, 75(14), 3631-3636.
 */

#ifndef __SABLIB_WHITTAKER_H__
#define __SABLIB_WHITTAKER_H__

#include <stdexcept>
#include <vector>

#include "../misc/diff.h"

namespace sablib {

/**
 * @brief Performs Whittaker smoothing (std::vector<double> version, with weights).
 *
 * @param y The input data to be smoothed.
 * @param w Weights for each data point.
 * @param lambda Smoothing parameter (larger values lead to more smoothing, but may flatten peaks).
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @return The smoothed data.
 * @throw std::invalid_argument One or more parameters wrong.
 */
const std::vector<double> Whittaker(
	const std::vector<double> & y, const std::vector<double> & w,
	const double lambda, const unsigned int s = 2
);

/**
 * @brief Performs Whittaker smoothing (std::vector<double> version, without weights).
 *
 * @param y The input data to be smoothed.
 * @param lambda Smoothing parameter (larger values lead to more smoothing, but may flatten peaks).
 * @param s The order of the difference (usually s = 1, 2, or 3).
 * @return The smoothed data.
 */
const std::vector<double> Whittaker(
	const std::vector<double> & y, const double lambda, const unsigned int s = 2
);

/**
 * @brief Performs Whittaker smoothing.
 *
 * @param y The input data to be smoothed.
 * @param w Weights for each data point.
 * @param lambdaDTD The matrix used for smoothness (lambda * D' * D).
 * @return The smoothed data.
 */
template <typename Derived>
const typename Derived::PlainObject
Whittaker(
	const Eigen::MatrixBase<Derived> & y,
	const Eigen::MatrixBase<Derived> & w,
	const Eigen::SparseMatrix<typename Derived::PlainObject::Scalar> & lambdaDTD
)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: y and w are not vector.");

	using Scalar = typename Derived::PlainObject::Scalar;

	Eigen::VectorXd z;
	Eigen::SparseMatrix<typename Derived::PlainObject::Scalar> W;
	Eigen::SimplicialCholesky< Eigen::SparseMatrix<Scalar> > solver;

	W = w.asDiagonal();
	solver.compute(W + lambdaDTD);

	if(solver.info() != Eigen::Success) {
		throw std::runtime_error("Whittaker(): solver calculation fails.");
	}

	z = solver.solve(W * y);

	return z;
}

}; // namespace sablib

#endif // __SABLIB_WHITTAKER_H__
