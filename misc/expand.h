/**
 * @file expand.h
 * @brief Expands the boundaries of the data
 * @author Izadori
 */

#ifndef __SABLIB_EXPAND_H__
#define __SABLIB_EXPAND_H__

#include <stdexcept>
#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief Expands the boundaries of a vector by padding with the first and last elements.
 *
 * @param y The input vector to be expanded.
 * @param n The number of elements to pad at each end.
 * @return The expanded vector with padded boundaries.
 */
template <typename Derived>
const typename Derived::PlainObject ExpandBoundaries(const Eigen::MatrixBase<Derived> & y, const unsigned int n)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: y is not vector.");

	typename Derived::PlainObject result(y.size() + 2 * n);

	if(y.cols() == 1) {
		result.block(0, 0, n, 1).fill(y(0));
		result.block(result.size() - n, 0, n, 1).fill(y(y.size() - 1));
		result.block(n, 0, y.size(), 1) = y;
	}
	else if(y.rows() == 1) {
		result.block(0, 0, 1, n).fill(y(0));
		result.block(0, result.size() - n, 1, n).fill(y(y.size() - 1));
		result.block(0, n, 1, y.size()) = y;
	}

	return result;
}

/**
 * @brief Trims the specified number of elements from both ends of a vector.
 *
 * @tparam Derived The derived type of the Eigen object.
 * @param y The input vector to be trimmed.
 * @param n The number of elements to remove from each end.
 * @return The trimmed vector.
 * @exception std::invalid_argument Thrown if n is larger than half the length of y.
 */
template <typename Derived>
const typename Derived::PlainObject TrimBoundaries(const Eigen::MatrixBase<Derived> & y, const unsigned int n)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived::IsVectorAtCompileTime, "Error: y is not vector.");

	if(y.size() < 2 * n) {
		throw std::invalid_argument("TrimBoubdaries(): n is too large.");
	}

	if(y.cols() == 1) {
		return y.block(n, 0, y.size() - 2 * n, 1);
	}
	else if(y.rows() == 1) {
		return y.block(0, n, 1, y.size() - 2 * n);
	}
}

}; // namespace sablib

#endif // __SABLIB_EXPAND_H__
