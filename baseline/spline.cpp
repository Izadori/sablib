/**
 * @file spline.cpp
 * @brief Baseline estimation from cubic spline(implementation)
 * @author Izadori
 */

#include "../misc/cubic_spline.h"

#include "spline.h"

namespace sablib {

//
// Implementation of BaselineSpline() function
//
const BaselineResult BaselineSpline(const std::vector<double> & y, const std::vector<unsigned int> & indices)
{
	if(y.size() == 0 || indices.size() == 0) {
		throw std::invalid_argument("BaselineSpline(): the length of y or indices is zero.");
	}

	if(y.size() < indices.size()) {
		throw std::invalid_argument("BaselineSpline(): the length of indices is larger than y.");
	}

	double max_index = y.size() - 1;
	Eigen::VectorXd xx(indices.size()), yy(indices.size());

	for(unsigned int i = 0; i < indices.size(); i++) {
		xx(i) = indices[i] / max_index;
		yy(i) = y[indices[i]];
	}

	CubicSpline spline(xx, yy);
	BaselineResult result;
	result.baseline.resize(y.size());
	result.corrected.resize(y.size());
	result.baseline = y;

	for(unsigned int i = indices[0]; i < indices.back(); i++) {
		result.baseline[i] = spline.Interpolate(i / max_index);
	}

	for(unsigned int i = 0; i < y.size(); i++) {
	    result.corrected[i] = y[i] - result.baseline[i];
	}

	return result;
}

}; // namespace sablib
