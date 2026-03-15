/**
 * @file spline.cpp
 * @brief Baseline estimation from cubic spline(implementation)
 * @author Izadori
 */

#include "spline.h"

namespace sablib {

//
// Implementation of BaselineSpline() function
//
const std::vector<double> BaselineSpline(const std::vector<double> & y, const std::vector<unsigned int> & indices)
{
	double max_index = y.size() - 1;
	Eigen::VectorXd xx(indices.size()), yy(indices.size());

	for(unsigned int i = 0; i < indices.size(); i++) {
		xx(i) = indices[i] / max_index;
		yy(i) = y[indices[i]];
	}

	CubicSpline spline(xx, yy);
	std::vector<double> result = y;

	for(unsigned int i = indices[0]; i < indices.back(); i++) {
		result[i] = spline.Interpolate(i / max_index);
	}

	return result;
}

}; // namespace sablib
