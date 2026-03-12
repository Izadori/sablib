//
// spdiags.h
//
// Copyright (c) 2026 Izadori
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef __SABLIB_SPDIAGS_H__
#define __SABLIB_SPDIAGS_H__

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief Returns a sparse matrix with the specified elements on its diagonals.
 * 
 * @param data  Matrix containing the elements to be placed on the diagonals as columns.
 * @param diags Vector indicating the positions where the elements should be placed (0: diagonal, positive: upper triangular, negative: lower triangular).
 * @param m     Number of rows in the resulting sparse matrix (if -1, determined by data).
 * @param n     Number of columns in the resulting sparse matrix (if -1, determined by data).
 * @return A sparse matrix with the elements placed on the specified diagonals.
 * @exception std::invalid_argument Thrown if the size of diags is larger than the number of rows in data.
 */
template <typename Derived>
Eigen::SparseMatrix<typename Derived::PlainObject::Scalar>
Spdiags(const Eigen::MatrixBase<Derived> & data, const Eigen::VectorXi & diags, const int m = -1, const int n = -1)
{
	using Scalar = typename Derived::PlainObject::Scalar;
	using T = Eigen::Triplet<Scalar>;

	if(diags.size() > data.rows()) {
		throw std::invalid_argument("diags size is larger than rows of data.");
	}

	int row_size, column_size;
	Eigen::SparseMatrix<Scalar> a;
	std::vector<T> triplets;

	row_size = (m <= 0) ? (int)data.rows() : m;
	column_size = (n <= 0) ? (int)data.cols() : n;

	a.resize(row_size, column_size);
	triplets.resize(row_size * diags.size());

	for (int k = 0; k < diags.size(); k++) {
		int start_index = std::max(0, diags(k));
		int end_index = std::min((int)data.cols(), column_size);

		for(int i = start_index; i < end_index; i++) {
			if(i - diags(k) < row_size && i < column_size) {
				triplets.push_back(T(i - diags(k), i, data(k, i)));
			}
		}
	}

	a.setFromTriplets(triplets.begin(), triplets.end());

	return a;
}

/**
 * @brief Returns a sparse matrix with the specified elements on its diagonals (std::vector<int> version).
 * 
 * @param data  Matrix containing the elements to be placed on the diagonals as columns.
 * @param diags Vector indicating the positions where the elements should be placed.
 * @param m     Number of rows in the resulting sparse matrix.
 * @param n     Number of columns in the resulting sparse matrix.
 * @return A sparse matrix with the elements placed on the specified diagonals.
 * @exception std::invalid_argument Thrown if the size of diags is larger than the number of rows in data.
 */
template <typename Derived>
inline Eigen::SparseMatrix<typename Derived::PlainObject::Scalar>
Spdiags(const Eigen::MatrixBase<Derived> & data, const std::vector<int> & diags, const int m = -1, const int n = -1)
{
	return Spdiags(
		data,
		Eigen::Map<Eigen::VectorXi>((int *)diags.data(), diags.size()),
		m, n
	);
}

/**
 * @brief Returns a sparse matrix with the specified elements on its diagonals (std::initializer_list<int> version).
 * 
 * @param data  Matrix containing the elements to be placed on the diagonals as columns.
 * @param diags List indicating the positions where the elements should be placed.
 * @param m     Number of rows in the resulting sparse matrix.
 * @param n     Number of columns in the resulting sparse matrix.
 * @return A sparse matrix with the elements placed on the specified diagonals.
 * @exception std::invalid_argument Thrown if the size of diags is larger than the number of rows in data.
 */
template <typename Derived>
inline Eigen::SparseMatrix<typename Derived::PlainObject::Scalar>
Spdiags(const Eigen::MatrixBase<Derived> & data, const std::initializer_list<int> & diags, const int m = -1, const int n = -1)
{
	return Spdiags(data, std::vector<int>{diags}, m, n);
}

}; // namespace sablib

#endif // __IZADORI_EIGEN_SPDIAGS_H__
