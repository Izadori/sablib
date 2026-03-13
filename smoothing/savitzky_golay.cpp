/**
 * @file savitzky_golay.cpp
 * @brief Smoothing using Savitzky-Golay filter(implementation)
 * @author Izadori
 */

#include <cmath>
#include <cassert>

#include "../misc/convolve.h"

#include "savitzky_golay.h"

namespace sablib {

//
// Implementation of SavitzkyGolayCoefficients() function
//
const std::vector<double> SavitzkyGolayCoefficients(
	const unsigned int n, const unsigned int polyorder, const unsigned derive, const double delta
)
{
	unsigned int window = 2 * n + 1;
	double p = 1;

	assert(window > polyorder && window > derive && polyorder > derive);

	Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(window, (int)-n, n);
	Eigen::MatrixXd x = Eigen::MatrixXd::Ones(window, polyorder + 1);

	for(unsigned int i = 1; i <= polyorder; i++){
		x.col(i) = (x.col(i - 1).array() * v.array()).matrix();
	}

	Eigen::MatrixXd coeff_mat = (x.transpose() * x).inverse() * x.transpose();

	if(derive > 0){
		for(unsigned int i = 1; i <= derive; i++){
			p *= i;
		}
		p /= std::pow(delta, derive);
	}

	Eigen::VectorXd coeff = coeff_mat.row(derive) * p;
	std::vector<double> result(coeff.size());

	Eigen::VectorXd::Map(result.data(), result.size()) = coeff;

	return result;
}

//
// Implementation of SavitzkyGolay() function
//
const std::vector<double> SavitzkyGolay(
	const std::vector<double> & y,
	const unsigned int n, const unsigned int polyorder, const unsigned derive, const double delta
)
{
	std::vector<double> coeff = SavitzkyGolayCoefficients(n, polyorder, derive, delta);
	Eigen::VectorXd w = Eigen::VectorXd::Map(coeff.data(), coeff.size());
	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());

	Eigen::VectorXd z = Convolve(yy, w, ConvolveMode::Same);

	std::vector<double> result(z.size());

	Eigen::VectorXd::Map(result.data(), result.size()) = z;

	return result;
}

}; // namespace sablib
