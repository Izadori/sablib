/**
 * @file bspline.h
 * @brief B-Spline basis functions and interpolation.
 * @author Izadori
 */

#ifndef __SABLIB_BSPLINE_H__
#define __SABLIB_BSPLINE_H__

#include <algorithm>
#include <type_traits>
#include <vector>
#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief A class for B-Spline operations including basis function calculation and interpolation.
 *
 * @tparam Scalar The scalar type for the spline (e.g., double, float).
 */
template <typename Scalar>
class BSpline final
{
	static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type.");

public:
	BSpline() = delete;

	/**
	 * @brief Constructor that initializes a B-Spline with a given degree and knot vector.
	 *
	 * @param degree The degree of the B-Spline.
	 * @param knots The knot vector.
	 */
	BSpline(const int degree, const Eigen::VectorX<Scalar> & knots);

	/**
	 * @brief Constructor that initializes a B-Spline and fits it to the given data points.
	 *
	 * @param degree The degree of the B-Spline.
	 * @param knots The knot vector.
	 * @param x The x-coordinates of the data points.
	 * @param y The y-coordinates of the data points.
	 */
	BSpline(
		const int degree, const Eigen::VectorX<Scalar> & knots,
		const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y
	);

	/**
	 * @brief Constructor that initializes a B-Spline with a given degree, knot vector, and coefficients.
	 *
	 * @param degree The degree of the B-Spline.
	 * @param knots The knot vector.
	 * @param coefficients The B-Spline coefficients.
	 */
	BSpline(
		const int degree, const Eigen::VectorX<Scalar> & knots, const Eigen::VectorX<Scalar> & coefficients
	);

	/**
	 * @brief Gets the number of basis functions.
	 *
	 * @return The number of basis functions.
	 */
	int BasisSize() const;

	/**
	 * @brief Returns the knot vector.
	 *
	 * @return The knot vector.
	 */
	const Eigen::VectorX<Scalar> Knots() const;

	/**
	 * @brief Returns the internal B-Spline coefficients.
	 *
	 * @return The B-Spline coefficients.
	 */
	const Eigen::VectorX<Scalar> Coefficients() const;

	/**
	 * @brief Constructs the design matrix (collocation matrix) for a given set of x-coordinates.
	 *
	 * @param x The x-coordinates of the data points.
	 * @return The sparse design matrix.
	 */
	const Eigen::SparseMatrix<Scalar> DesignMatrix(const Eigen::VectorX<Scalar> & x) const;

	/**
	 * @brief Fits the B-Spline to the given data points by calculating the coefficients.
	 *
	 * @param x The x-coordinates of the data points.
	 * @param y The y-coordinates of the data points.
	 */
	void Fit(const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y);

	/**
	 * @brief Interpolates the value at a given x-coordinate using provided B-Spline coefficients.
	 *
	 * @param x The x-coordinate to interpolate.
	 * @param coefficients The B-Spline coefficients.
	 * @return The interpolated value at x.
	 */
	const Scalar Interpolate(const Scalar x, const Eigen::VectorX<Scalar> & coefficients) const;

	/**
	 * @brief Interpolates the value at a given x-coordinate using the internal coefficients.
	 *
	 * @param x The x-coordinate to interpolate.
	 * @return The interpolated value at x.
	 */
	const Scalar Interpolate(const Scalar x) const;

private:
	int sp_degree;                          /**< Degree of the B-Spline. */
	Eigen::VectorX<Scalar> sp_knots;        /**< Knot vector of the B-Spline. */
	int basis;                              /**< Number of basis functions. */
	Eigen::VectorX<Scalar> sp_coefficients; /**< Coefficients (control points) of the B-Spline. */

	/**
	 * @brief Finds the knot span for a given value x.
	 *
	 * @param x The value to find the span for.
	 * @return The index of the knot span.
	 */
	int FindSpan(const Scalar x) const;

	/**
	 * @brief Calculates the non-zero basis functions at a given x-coordinate and knot span.
	 *
	 * @param span The index of the knot span.
	 * @param x The x-coordinate to evaluate.
	 * @return A vector of non-zero basis function values.
	 */
	const Eigen::VectorX<Scalar> BasisFunctions(const int span, const Scalar x) const;
};

//
// Implementation of constructor
//
template <typename Scalar>
BSpline<Scalar>::BSpline(
	const int degree, const Eigen::VectorX<Scalar> & knots
) : sp_degree(degree), sp_knots(knots), basis(knots.size() - degree - 1)
{
}

//
// Implementation of constructor with fitting
//
template <typename Scalar>
BSpline<Scalar>::BSpline(
	const int degree, const Eigen::VectorX<Scalar> & knots,
	const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y
) : sp_degree(degree), sp_knots(knots), basis(knots.size() - degree - 1)
{
	Fit(x, y);
}

//
// Implementation of constructor with coefficients
//
template <typename Scalar>
BSpline<Scalar>::BSpline(
	const int degree, const Eigen::VectorX<Scalar> & knots, const Eigen::VectorX<Scalar> & coefficients
) : sp_degree(degree), sp_knots(knots), basis(knots.size() - degree - 1), sp_coefficients(coefficients)
{
}

//
// Implementation of BasisSize() method
//
template <typename Scalar>
int BSpline<Scalar>::BasisSize() const
{
	return basis;
}

//
// Implementation of Knots() method
//
template <typename Scalar>
const Eigen::VectorX<Scalar> BSpline<Scalar>::Knots() const
{
	return sp_knots;
}

//
// Implementation of Coefficients() method
//
template <typename Scalar>
const Eigen::VectorX<Scalar> BSpline<Scalar>::Coefficients() const
{
	return sp_coefficients;
}

//
// Implementation of FindSpan() method
//
template <typename Scalar>
int BSpline<Scalar>::FindSpan(const Scalar x) const
{
	int n = basis - 1;

	if (x >= sp_knots(n + 1)) {
		return n;
	}

	if (x <= sp_knots(sp_degree)) {
		return sp_degree;
	}

	int low  = sp_degree;
	int high = n + 1;
	int mid  = (low + high) / 2;

	while (x < sp_knots(mid) || x >= sp_knots(mid + 1)) {
		if (x < sp_knots(mid)) {
			high = mid;
		}
		else {
			low = mid;
		}

		mid = (low + high) / 2;
	}

	return mid;
}

//
// Implementation of BasisFunctions() method
//
template <typename Scalar>
const Eigen::VectorX<Scalar> BSpline<Scalar>::BasisFunctions(const int span, const Scalar x) const
{
	Eigen::VectorX<Scalar> N(sp_degree + 1);
	Eigen::VectorX<Scalar> left(sp_degree + 1);
	Eigen::VectorX<Scalar> right(sp_degree + 1);

	N(0) = 1.0;

	for(int j = 1; j <= sp_degree; j++) {
		left(j)  = x - sp_knots(span + 1 - j);
		right(j) = sp_knots(span + j) - x;

		Scalar saved = 0.0;

		for(int r = 0; r < j; r++) {
			Scalar temp = N(r) / (right(r + 1) + left(j - r));
			N(r) = saved + right(r + 1) * temp;
			saved = left(j - r) * temp;
		}

		N(j) = saved;
	}

	return N;
}

//
// Implementation of DesignMatrix() method
//
template <typename Scalar>
const Eigen::SparseMatrix<Scalar> BSpline<Scalar>::DesignMatrix(const Eigen::VectorX<Scalar> & x) const
{
	int n = x.size();
	std::vector< Eigen::Triplet<Scalar> > triplets;

	for(int r = 0; r < n; r++) {
		Scalar xi = x(r);
		int span = FindSpan(xi);
		Eigen::VectorX<Scalar> N = BasisFunctions(span, xi);

		for(int j = 0; j <= sp_degree; j++) {
			int col = span - sp_degree + j;

			if(col >=0 && col < basis) {
				triplets.emplace_back(r, col, N(j));
			}
		}
	}

	Eigen::SparseMatrix<Scalar> B(n , basis);
	B.setFromTriplets(triplets.begin(), triplets.end());

	return B;
}

//
// Implementation of Fit() method
//
template <typename Scalar>
void BSpline<Scalar>::Fit(const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y)
{
	auto B = DesignMatrix(x);
	Eigen::SparseMatrix<Scalar> BTB = B.transpose() * B;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix<Scalar> > solver;

	solver.compute(BTB);

	if(solver.info() != Eigen::Success) {
		throw std::runtime_error("BSpline::Fit(): solver calculation fails.");
	}

	sp_coefficients = solver.solve(B.transpose() * y);
}

//
// Implementation of Interpolate() method
//
template<typename Scalar>
const Scalar BSpline<Scalar>::Interpolate(const Scalar x, const Eigen::VectorX<Scalar> & coefficients) const
{
	int span = FindSpan(x);
	Eigen::VectorX<Scalar> N = BasisFunctions(span, x);

	Scalar y = 0.0;

	for(int j = 0; j <= sp_degree; j++) {
		int idx = span - sp_degree + j;

		if(idx >=0 && idx < coefficients.size()) {
			y += N(j) * coefficients(idx);
		}
	}

	return y;
}

//
// Implementation of Interpolate() method with internal coefficients
//
template<typename Scalar>
inline const Scalar BSpline<Scalar>::Interpolate(const Scalar x) const
{
	return Interpolate(x, sp_coefficients);
}

} // namespace sablib

#endif // __SABLIB_BSPLINE_H__
