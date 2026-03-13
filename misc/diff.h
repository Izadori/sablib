/**
 * @file diff.h
 * @brief MATLAB-like diff() function
 * @author Izadori
 */

#ifndef __SABLIB_DIFF_H__
#define __SABLIB_DIFF_H__

#include <Eigen/Eigen>

namespace sablib {

/**
 * @brief Direction for difference calculation.
 */
enum class Dir {
	RowWise,   /**< Calculates the difference along the rows (MATLAB-like).*/
	ColumnWise /**< Calculates the difference along the columns (Python-like).*/
};

/**
 * @brief Calculates the n-th discrete difference along the given axis.
 *
 * @param m0  The input matrix or vector.
 * @param n   The order of difference (default is 1).
 * @param dir The direction to calculate the difference (RowWise or ColumnWise).
 * @return The n-th discrete difference.
 */
template <typename Derived>
const typename Derived::PlainObject
Diff(const Eigen::MatrixBase<Derived> & m0, const int n = 1, const Dir dir = Dir::RowWise)
{
	typename Derived::PlainObject m = m0;
	int rows = m.rows(), columns = m.cols();

	if(dir == Dir::ColumnWise){
		for(int i = 0; i < n ; i++){
			int col_counts = columns - i - 1;
			m.leftCols(col_counts) = m.block(0, 1, rows, col_counts) - m.leftCols(col_counts);
		}

		return m.leftCols(columns - n);
	}
	else{
		for(int i = 0; i < n ; i++){
			int row_counts = rows - i - 1;
			m.topRows(row_counts) = m.block(1, 0, row_counts, columns) - m.topRows(row_counts);
		}

		return m.topRows(rows - n);
	}
}

/**
 * @brief Calculates the n-th discrete difference for sparse matrices along the given axis.
 *
 * @param m0  The input sparse matrix or vector.
 * @param n   The order of difference (default is 1).
 * @param dir The direction to calculate the difference (RowWise or ColumnWise).
 * @return The n-th discrete difference.
 */
template <typename Derived>
const typename Derived::PlainObject
Diff(const Eigen::SparseMatrixBase<Derived> & m0, const int n = 1, const Dir dir = Dir::RowWise)
{
	typename Derived::PlainObject m = m0;

	if(dir == Dir::ColumnWise){
		for(int i = 0; i < n; i++){
			m = (m.rightCols(m.cols() - 1) - m.leftCols(m.cols() - 1)).eval();
		}
	}
	else{
		for(int i = 0; i < n; i++){
			m = (m.bottomRows(m.rows() - 1) - m.topRows(m.rows() - 1)).eval();
		}
	}

	return m;
}

}; // namespace sablib

#endif // __SABLIB_DIFF_H__
