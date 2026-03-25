/**
 * @file polyfit.h
 * @brief Polynomial fitting using the least squares method (Gauss-Newton for linear models) and evaluating polynomial function
 * @author Izadori
 */

#ifndef __SABLIB_POLYFIT_H__
#define __SABLIB_POLYFIT_H__

#include <stdexcept>
#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief Generates a Vandermonde matrix for a given vector and polynomial order.
 *
 * @param x The input vector of data points.
 * @param polyorder The order of the polynomial (columns in the matrix will be polyorder + 1).
 * @return The generated Vandermonde matrix.
 * @exception std::invalid_argument If the length of x is zero or polyorder is zero.
 */
template <typename Derived>
const Eigen::MatrixX<typename Derived::PlainObject::Scalar>
Vandermonde(const Eigen::MatrixBase<Derived> & x, const unsigned int polyorder)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: x is not vector.");

	if(x.size() == 0) {
		throw std::invalid_argument("Vandermonde(): the length of x is zero.");
	}

	if(polyorder == 0) {
		throw std::invalid_argument("Vandermonde(): polyorder is zero.");
	}

	using Scalar = typename Derived::PlainObject::Scalar;
	unsigned int n = x.size();
	unsigned int m = polyorder + 1;
	Eigen::MatrixX<Scalar> V(n, m);

	for(int i = 0; i < n; ++i) {
		Scalar val = 1.0;

		for(int j = 0; j < m; ++j) {
			V(i, j) = val;
			val *= x(i);
		}
	}

	return V;
}

/**
 * @brief Solves the polynomial least squares fitting problem using a pre-calculated Vandermonde matrix.
 *
 * @param V The Vandermonde matrix.
 * @param y The target vector of data points.
 * @return The coefficients of the fitted polynomial.
 * @exception std::invalid_argument If the number of rows in V does not match the size of y.
 */
template <typename Derived>
const typename Derived::PlainObject
PolyFit(
	const Eigen::MatrixX<typename Derived::PlainObject::Scalar> & V,
	const Eigen::MatrixBase<Derived> & y
)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: y is not vector.");

	if(V.rows() != y.size()) {
		throw std::invalid_argument("PolyFit(): the sizes of y and V do not match.");
	}

	return (V.transpose() * V).ldlt().solve(V.transpose() * y);
}

/**
 * @brief Fits a polynomial of a specified order to the given data points.
 *
 * This function uses the normal equations, which is the Gauss-Newton method applied to a linear model.
 *
 * @param x The x-coordinates of the data points.
 * @param y The y-coordinates of the data points.
 * @param polyorder The order of the polynomial to fit.
 * @return The coefficients of the fitted polynomial (from lowest to highest order).
 * @exception std::invalid_argument If x and y sizes differ or if polyorder is too high for the number of points.
 */
template <typename Derived>
const typename Derived::PlainObject
PolyFit(const Eigen::MatrixBase<Derived> & x, const Eigen::MatrixBase<Derived> & y, const unsigned int polyorder)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: x or y is not vector.");

	if(x.size() != y.size()) {
		throw std::invalid_argument("PolyFit(): x and y must have the same size.");
	}

	if(polyorder + 1 > x.size()) {
		throw std::invalid_argument("PolyFit(): polyorder is too high for the number of points.");
	}

	return PolyFit(Vandermonde(x, polyorder), y);
}

/**
 * @brief Evaluates the polynomial at specified points using a pre-calculated Vandermonde matrix.
 *
 * @param coeff The coefficients of the polynomial.
 * @param V The Vandermonde matrix.
 * @return The evaluated values.
 * @exception std::invalid_argument If the length of coeff is zero or the number of columns in V does not match the size of coeff.
 */
template <typename Derived>
const typename Derived::PlainObject
PolyVal(
	const Eigen::MatrixBase<Derived> & coeff,
	const Eigen::MatrixX<typename Derived::PlainObject::Scalar> & V
)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: coeff is not vector.");

	if(coeff.size() == 0) {
		throw std::invalid_argument("PolyVal(): the length of coeff is zero.");
	}

	if(V.cols() != coeff.size()) {
		throw std::invalid_argument("PolyVal(): the sizes of V and coeff do not match.");
	}

	return V * coeff;
}

/**
 * @brief Evaluates the polynomial at specified x-coordinates.
 *
 * @param coeff The coefficients of the polynomial.
 * @param x The x-coordinates where the polynomial will be evaluated.
 * @return The evaluated values.
 */
template <typename Derived>
const typename Derived::PlainObject
PolyVal(const Eigen::MatrixBase<Derived> & coeff, const Eigen::MatrixBase<Derived> & x)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: coeff or x is not vector.");

	return PolyVal(coeff, Vandermonde(x, coeff.size() - 1));
}

/**
 * @brief Evaluates the polynomial at single specified x-coordinate.
 *
 * @param coeff The coefficients of the polynomial.
 * @param x The x-coordinate where the polynomial will be evaluated.
 * @return The evaluated value.
 * @exception std::invalid_argument If the length of coeff is zero.
 */
template <typename Derived>
const typename Derived::PlainObject::Scalar
PolyVal(const Eigen::MatrixBase<Derived> & coeff, const typename Derived::PlainObject::Scalar x)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: coeff is not vector.");

	if(coeff.size() == 0) {
		throw std::invalid_argument("PolyVal(): the length of coeff is zero.");
	}

	typename Derived::PlainObject::Scalar val = 0, xx = 1.0;

	for(int i = 0; i < coeff.size(); i++) {
		val += coeff(i) * xx;
		xx  *= x;
	}

	return val;
}

}; // namespace sablib

#endif // __SABLIB_POLYFIT_H__
