/**
 * @file rubberband.cpp
 * @brief Baseline estimation using Rubber Band method(implementation)
 * @author Izadori
 */

#include <Eigen/Dense>

#include "rubberband.h"

namespace sablib {

//
// Implementation of BaselineRubberBand() function
//
const BaselineResult BaselineRubberBand(const std::vector<double> & y)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineRubberBand(): the length of y is zero.");
	}

	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());
	Eigen::VectorXd z = Eigen::VectorXd::Zero(yy.size());

	// Lower convex hull
	std::vector<unsigned int> lower_hull;

	lower_hull.reserve(yy.size());

	for(unsigned int i = 0; i < yy.size(); i++) {
		while(lower_hull.size() >= 2) {
			unsigned int j = lower_hull[lower_hull.size() - 2];
			unsigned int k = lower_hull[lower_hull.size() - 1];

			if((k - j) * (yy(i) - yy(j)) - (i - j) * (yy(k) - yy(j)) <= 0) {
				lower_hull.pop_back();
			}
			else {
				break;
			}
		}

		lower_hull.push_back(i);
	}

    for (size_t h = 0; h < lower_hull.size() - 1; ++h)
    {
        int i0 = lower_hull[h];
        int i1 = lower_hull[h + 1];

        double y0 = yy(i0);
        double y1 = yy(i1);

        for (int i = i0; i <= i1; ++i)
        {
            double t = double(i - i0) / (i1 - i0);
            z(i) = y0 + t * (y1 - y0);
        }
    }

	BaselineResult result;
	result.baseline.resize(z.size());
	result.corrected.resize(z.size());

	Eigen::VectorXd::Map(result.baseline.data(), z.size()) = z;
	Eigen::VectorXd::Map(result.corrected.data(), z.size()) = yy - z;

	return result;
}

}; // namespace sablib
