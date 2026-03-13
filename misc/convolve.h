/**
 * @file convolve.h
 * @brief Python's NumPy-like convolve() function
 * @author Izadori
 */

#ifndef __SABLIB_CONVOLVE_H__
#define __SABLIB_CONVOLVE_H__

#include <algorithm>
#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief Mode of Convolve() function
 */
enum class ConvolveMode
{
	Full, /**< This returns the convolution at each point of overlap, with an output shape of `N + M - 1`.*/
	Same, /**< This returns output of length `max(M, N)`.*/
	Valid /**< This returns output of length `max(M, N) - min(M, N) + 1`. The convolution product is only given for points where the signals overlap completely. */
};

/**
 * @brief Returns the discrete, linear convolution of two one-dimensional sequences.
 *
 * @param v1 First one-dimensional input array.
 * @param v2 Second one-dimensional input array.
 * @param mode The mode of the convolution (default is Full).
 * @return The convolved one-dimensional sequence.
 */
template <typename Derived1, typename Derived2>
const typename Derived1::PlainObject
Convolve(
	const Eigen::MatrixBase<Derived1> & v1,
	const Eigen::MatrixBase<Derived2> & v2,
	const ConvolveMode mode = ConvolveMode::Full
)
{
	// Although parameters are received as MatrixBase<Derived>, only vector classes are allowed.
	// Others will be rejected at compile time.
	static_assert(Derived1::IsVectorAtCompileTime, "Error: v1 is not vector.");
	static_assert(Derived2::IsVectorAtCompileTime, "Error: v2 is not vector.");

	int length, start_index, max_size, min_size;
	typename Derived1::PlainObject result;

	max_size = std::max(v1.size(), v2.size());
	min_size = std::min(v1.size(), v2.size());

	length = v1.size() + v2.size() - 1;
	start_index = 0;
	result = Derived1::PlainObject::Zero(length);

	for(int j = 0; j < v2.size(); j++) {
		for(int i = 0; i < v1.size(); i++) {
			result(i + j) += v1(i) * v2(j);
		}
	}

	if (mode == ConvolveMode::Valid) {
		length = max_size - min_size + 1;
		start_index = min_size - 1;
	}
	else if(mode == ConvolveMode::Same) {
		length = max_size;
		start_index = (int)(min_size / 2.0 + 0.5) - 1;
	}

	if(result.rows() == 1) {
		return result.block(0, start_index, 1, length);
	}
	else {
		return result.block(start_index, 0, length, 1);
	}
}

}; // namespace sablib

#endif // __SABLIB_CONVOLVE_H__
