/**
 * @file polynomial.cpp
 * @brief Baseline estimation with polynomial line(implementation)
 * @author Izadori
 */

#include "polynomial.h"

namespace sablib {

//
// Implementation of BaselineLinear() function
//
std::vector<double> BaselineLinear(std::vector<double> & y, const unsigned int index1, const unsigned int index2)
{
	if(index1 >= index2 || y.size() <= index1 || y.size() <= index2) {
		throw std::invalid_argument("BaselineLinear(): illegal indices");
	}

	std::vector<double> result = y;
	double m = (y[index2] - y[index1]) / (index2 - index1);

	auto f = [&](const unsigned int x) {
		return m * (x - index1) + y[index1];
	};

	for(unsigned int x = index1; x <= index2; x++) {
		result[x] = f(x);
	}

	return result;
}

//
// Implementation of BaselinePolynomial() function
//
std::vector<double> BaselinePolynomial(
	std::vector<double> & y, const unsigned int polyorder, const std::vector<unsigned int> & indices
)
{
	if(polyorder >= indices.size()) {
		throw std::invalid_argument("BaselinePolynomial(): Too few indices.");
	}

	double max_index = y.size();
	std::vector<unsigned int> sorted_indices = indices;
	Eigen::VectorXd xx(indices.size()), yy(indices.size());

	std::sort(sorted_indices.begin(), sorted_indices.end());

	for(unsigned int i = 0; i < indices.size(); i++) {
		xx(i) = sorted_indices[i] / max_index;
		yy(i) = y[sorted_indices[i]];
	}

	Eigen::VectorXd coefficients = PolyFit(xx, yy, polyorder);

	std::vector<double> result = y;

	for(unsigned int i = sorted_indices[0]; i < sorted_indices.back(); i++) {
		double x = i / max_index;
		double fx = 0;
		double k = 1;

		for(int j = 0; j < coefficients.size(); j++) {
			fx += coefficients(j) * k;
			k *= x;
		}

		result[i] = fx;
	}

	return result;
}

}; // namespace sablib
