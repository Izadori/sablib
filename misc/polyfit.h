/**
 * @file polyfit.h
 * @brief Polynomial fitting using the least squares method (Gauss-Newton for linear models)
 * @author Izadori
 */

#ifndef __SABLIB_POLYFIT_H__
#define __SABLIB_POLYFIT_H__

#include <stdexcept>
#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief Fits a polynomial of a specified order to the given data points.
 *
 * This function uses the normal equations, which is the Gauss-Newton method applied to a linear model.
 *
 * @tparam Derived The derived type of the Eigen object.
 * @param x The x-coordinates of the data points.
 * @param y The y-coordinates of the data points.
 * @param polyorder The order of the polynomial to fit.
 * @return The coefficients of the fitted polynomial (from lowest to highest order).
 * @exception std::invalid_argument Thrown if x and y sizes differ or if polyorder is too high.
 */
template <typename Derived>
const typename Derived::PlainObject
PolyFit(const Eigen::MatrixBase<Derived> & x, const Eigen::MatrixBase<Derived> & y, const unsigned int polyorder)
{
    using Scalar = typename Derived::PlainObject::Scalar;
    const int n = x.size();
    const int m = polyorder + 1;

    if (x.size() != y.size()) {
        throw std::invalid_argument("PolyFit(): x and y must have the same size.");
    }

    if (m > n) {
        throw std::invalid_argument("PolyFit(): polyorder is too high for the number of points.");
    }

    Eigen::MatrixX<Scalar> A(n, m);

	for (int i = 0; i < n; ++i) {
        Scalar val = 1.0;
        for (int j = 0; j < m; ++j) {
            A(i, j) = val;
            val *= x(i);
        }
    }

    return (A.transpose() * A).ldlt().solve(A.transpose() * y);
}

}; // namespace sablib

#endif // __SABLIB_POLYFIT_H__
