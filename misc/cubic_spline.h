/**
 * @file cubic_spline.h
 * @brief Performs cubic spline parameter calculation and interpolation.
 * @author Izadori
 */

#ifndef __SABLIB_CUBIC_SPLINE_H__
#define __SABLIB_CUBIC_SPLINE_H__

#include <stdexcept>
#include <vector>
#include <Eigen/Sparse>

namespace sablib {

/**
 * @brief A class for cubic spline interpolation.
 *
 * @tparam Scalar The scalar type for the spline (e.g., double, float).
 */
template <typename Scalar>
class CubicSpline final
{
public:
	/**
	 * @brief Default constructor.
	 */
	CubicSpline() = default;

	/**
	 * @brief Constructor that fits a cubic spline to the given data points.
	 *
	 * @param x The x-coordinates of the data points.
	 * @param y The y-coordinates of the data points.
	 */
	CubicSpline(const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y);

	/**
	 * @brief Fits a cubic spline to the given data points.
	 *
	 * @param x The x-coordinates of the data points.
	 * @param y The y-coordinates of the data points.
	 */
	void Fit(const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y);

	/**
	 * @brief Interpolates the value at a given x-coordinate.
	 *
	 * @param x The x-coordinate to interpolate.
	 * @return The interpolated value at x.
	 */
	Scalar Interpolate(const double x);

private:
	/**
	 * @brief Solves a tridiagonal system of equations using Thomas' algorithm (TDMA).
	 *
	 * @param lower The lower diagonal (size n). lower[0] is ignored.
	 * @param diag The main diagonal (size n).
	 * @param upper The upper diagonal (size n). upper[n-1] is ignored.
	 * @param rhs The right-hand side vector (size n).
	 * @return The solution vector.
	 */
	Eigen::VectorX<Scalar> SolveTridiagonal(
		const Eigen::VectorX<Scalar>& lower, const Eigen::VectorX<Scalar>& diag,
		const Eigen::VectorX<Scalar>& upper, const Eigen::VectorX<Scalar>& rhs
	);

	Eigen::VectorX<Scalar> a, b, c, d;
	Eigen::VectorX<Scalar> sp_x, sp_y;
};

//
// Implementation of constructor
//
template <typename Scalar>
inline CubicSpline<Scalar>::CubicSpline(const Eigen::VectorX<Scalar>& x, const Eigen::VectorX<Scalar>& y)
{
	Fit(x, y);
}

//
// Implementation of Fit() method
//
template <typename Scalar>
void CubicSpline<Scalar>::Fit(const Eigen::VectorX<Scalar> & x, const Eigen::VectorX<Scalar> & y)
{
	sp_x = x;
	sp_y = y;

	int n = x.size() - 1;
	Eigen::VectorX<Scalar> h(n);

	for(int i = 0; i < n; i++) {
		h(i) = x(i + 1) - x(i);
	}

	Eigen::VectorX<Scalar> lower = Eigen::VectorX<Scalar>::Zero(n + 1);
	Eigen::VectorX<Scalar> diag  = Eigen::VectorX<Scalar>::Ones(n + 1);
	Eigen::VectorX<Scalar> upper = Eigen::VectorX<Scalar>::Zero(n + 1);
	Eigen::VectorX<Scalar> rhs   = Eigen::VectorX<Scalar>::Zero(n + 1);

	for(int i = 1; i < n; i++) {
		lower(i) = h(i - 1);
		diag(i)  = 2 * (h(i - 1) + h(i));
		upper(i) = h(i);
		rhs(i)   = 6 * ((y(i + 1) - y(i)) / h(i) - (y(i) - y(i - 1)) / h(i - 1));
	}

	Eigen::VectorX<Scalar> m = SolveTridiagonal(lower, diag, upper, rhs);

	a.resize(n);
	b.resize(n);
	c.resize(n);
	d.resize(n);

	for(int i=0;i<n;i++) {
		a(i) = y(i);
		b(i) = (y(i + 1) - y(i)) / h(i) - h(i) * (2 * m(i) + m(i + 1)) / 6.0;
		c(i) = m(i) / 2.0;
		d(i) = (m(i + 1) - m(i)) / (6.0 * h(i));
	}
}

//
// Implementation of Thomas' algorithm (TDMA)
//
template <typename Scalar>
Eigen::VectorX<Scalar> CubicSpline<Scalar>::SolveTridiagonal(
	const Eigen::VectorX<Scalar>& lower, const Eigen::VectorX<Scalar>& diag,
	const Eigen::VectorX<Scalar>& upper, const Eigen::VectorX<Scalar>& rhs
)
{
	int n = diag.size();
	Eigen::VectorX<Scalar> cp = Eigen::VectorX<Scalar>::Zero(n);
	Eigen::VectorX<Scalar> dp = Eigen::VectorX<Scalar>::Zero(n);
	Eigen::VectorX<Scalar> x  = Eigen::VectorX<Scalar>::Zero(n);

	// Forward elimination
	cp(0) = upper(0) / diag(0);
	dp(0) = rhs(0) / diag(0);

	for (int i = 1; i < n; i++) {
		Scalar m = diag(i) - lower(i) * cp(i - 1);
		cp(i) = upper(i) / m;
		dp(i) = (rhs(i) - lower(i) * dp(i - 1)) / m;
	}

	// Backward substitution
	x(n - 1) = dp(n - 1);
	for (int i = n - 2; i >= 0; i--) {
		x(i) = dp(i) - cp(i) * x(i + 1);
	}

	return x;
}

//
// Implementation of Interpolate() method
//
template <typename Scalar>
Scalar CubicSpline<Scalar>::Interpolate(const double x)
{
	int n = sp_x.size() - 1;
	int i = n - 1;

	for(int j = 0; j < n; j++) {
		if(x >= sp_x(j) && x <= sp_x(j + 1)) {
			i = j;
			break;
		}
	}

	double dx = x - sp_x(i);

	return a(i) + b(i) * dx + c(i) * dx * dx + d(i) * dx * dx * dx;
}

}; // namespace sablib

#endif // __SABLIB_CUBIC_SPLINE_H__
