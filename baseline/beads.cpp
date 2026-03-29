/**
 * @file beads.cpp
 * @brief Baseline estimation and subtraction using Baseline Estimation And Denoising using Sparsity(BEADS)(implementation)
 * @author Izadori
 * @note This file includes code derived from [Duval's MATLAB implemetation](https://jp.mathworks.com/matlabcentral/fileexchange/49974-beads-baseline-estimation-and-denoising-with-sparsity/files/license.txt).
 *       Original MATLAB code is licensed under BSD 3-Clause License.
 */

/*
Original License:
Copyright (c) 2018, Laurent Duval
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution
* Neither the name of IFP Energies nouvelles nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <functional>

#include "../misc/convolve.h"
#include "../misc/diff.h"
#include "../misc/spdiags.h"

#include "beads.h"

namespace sablib {

namespace {

const std::tuple< Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double> >
BAfilt(const unsigned int d, const double frequency, const unsigned int length)
{
	Eigen::VectorXd a(1), b, b1(2), v(3), v2(2);
	double omega_c = 2 * M_PI * frequency;
	double t;

	b1 << 1, -1;
	v  << -1, 2, -1;

	for(unsigned int i = 1; i < d; i++){
		b1 = Convolve(b1, v).eval();
	}

	v2 << -1, 1;
	b = Convolve(b1, v2);

	v << 1, 2, 1;
	a << 1;

	for(unsigned int i = 0; i < d; i++){
		a = Convolve(a, v).eval();
	}

	t = std::pow((1 - std::cos(omega_c)) / (1 + std::cos(omega_c)), d);
	a = (b + t * a).eval();

	Eigen::MatrixXd xa(a.size(), length);
	Eigen::MatrixXd xb(b.size(), length);

	for(unsigned int i = 0; i < length; i++) {
		xa.block(0, i, a.size(), 1) = a;
		xb.block(0, i, b.size(), 1) = b;
	}

	Eigen::VectorXi dr = Eigen::VectorXi::LinSpaced(2 * d + 1, -d, d);
	Eigen::SparseMatrix<double> A = Spdiags(xa, dr, length, length);
	Eigen::SparseMatrix<double> B = Spdiags(xb, dr, length, length);

	return std::tuple< Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double> >(A, B);
}

}; // unnamed namespace

//
// Implementation of BaselineBeads() function
//
const std::tuple< std::vector<double>, std::vector<double> >
BaselineBeads(
	const std::vector<double> & y, const unsigned int s, const double frequency, const double r,
	const double lambda0, const double lambda1, const double lambda2,
	const unsigned int loop,  const double eps, const BeadsPenalty penalty
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineBeads(): the length of y is zero.");
	}

	if(s == 0 || s > 3) {
		throw std::invalid_argument("BaselineBeads(): s must be 1, 2 or 3.");
	}

	if(frequency <= 0) {
		throw std::invalid_argument("BaselineBeads(): non-positive frequency value is given.");
	}

	if(r <= 0) {
		throw std::invalid_argument("BaselineBeads(): non-positive r value is given.");
	}

	if(lambda0 <= 0 || lambda1 <= 0 || lambda2 <= 0) {
		throw std::invalid_argument("BaselineBeads(): non-positive lambda value is given.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineBeads(): loop is zero.");
	}

	if(eps <= 0) {
		throw std::invalid_argument("BaselineBeads(): non-positive eps value is given.");
	}

	const double eps0 = 1e-6;
	const double eps1 = 1e-6;

	Eigen::VectorXd yy = Eigen::VectorXd::Map(y.data(), y.size());

	std::function<double(const double)> phi, wfun;

	switch(penalty) {
		case BeadsPenalty::L1_v1:
			phi = [&](const double xx) {
					double abs_x = std::fabs(xx);
					return std::sqrt(abs_x * abs_x + eps1);
			};
			wfun = [&](const double xx) {
					return 1 / phi(xx);
			};
			break;
		case BeadsPenalty::L1_v2:
			phi = [&](const double xx) {
				double abs_x = std::fabs(xx);
				return abs_x - eps1 * std::log(abs_x + eps1);
			};
			wfun = [&](const double xx) {
				return 1 / (std::fabs(xx) + eps1);
			};
			break;
		default:
			throw std::invalid_argument("invalid penalty function type.");
	}

	auto theta = [&](const double xx) {
		if(xx > eps0) {
			return xx;
		}
		else if(xx < -eps0) {
			return -r * xx;
		}
		else {
			return (1 + r) * xx * xx / (4 * eps0) + (1 - r) * xx / 2 + eps0 * (1 + r) / 4;
		}
	};

	int length = yy.size();
	auto [ A, B ] = BAfilt(s, frequency, length);

	Eigen::SparseMatrix<double> I, D1, D2;

	I.resize(length, length);
	I.setIdentity();
	D1 = Diff(I);
	D2 = Diff(I, 2);

	Eigen::SparseLU< Eigen::SparseMatrix<double> > solverA, solverQ;

	solverA.compute(A);

	if(solverA.info() != Eigen::Success) {
		throw std::runtime_error("BaselineBeads(): solverA calculation fails.");
	}

	Eigen::SparseMatrix<double> BTB = B.transpose() * B;
	Eigen::VectorXd b = Eigen::VectorXd::Constant(length, (1 - r) / 2);
	Eigen::VectorXd d = BTB * (solverA.solve(yy)) - lambda0 * A.transpose() * b;
	Eigen::VectorXd w1 = Eigen::VectorXd::Constant(length - 1, lambda1);
	Eigen::VectorXd w2 = Eigen::VectorXd::Constant(length - 2, lambda2);

	Eigen::VectorXd x = yy;

	double prev_c = 0.5
		* (B * solverA.solve(x)).array().square().sum()
		+ lambda0 * x.unaryExpr(theta).sum()
		+ lambda1 * (D1 * x).unaryExpr(phi).sum()
		+ lambda2 * (D2 * x).unaryExpr(phi).sum();

	for(unsigned int i = 0; i < loop; i++) {
		Eigen::SparseMatrix<double> L1, L2, G, M;

		L1 = (w1.array() * (D1 * x).unaryExpr(wfun).array()).matrix().asDiagonal();
		L2 = (w2.array() * (D2 * x).unaryExpr(wfun).array()).matrix().asDiagonal();

		G = x.unaryExpr([&](const double xx) {
			double z = (-eps0 <= xx && xx <= eps0) ? eps0 : std::fabs(xx);
			return (1 + r) / 4 / z;
		}).asDiagonal();

		M = 2 * lambda0 * G + D1.transpose() * L1 * D1 + D2.transpose() * L2 * D2;
		solverQ.compute(BTB + A.transpose() * M * A);

		if(solverQ.info() != Eigen::Success) {
			throw std::runtime_error("BaselineBeads(): solverQ calculation fails.");
		}

		x = A * solverQ.solve(d);

		Eigen::VectorXd a = yy - x;
		double c = 0.5
			* (B * solverA.solve(a)).array().square().sum()
			+ lambda0 * x.unaryExpr(theta).sum()
			+ lambda1 * (D1 * x).unaryExpr(phi).sum()
			+ lambda2 * (D2 * x).unaryExpr(phi).sum();

		if(std::fabs((prev_c - c) / c) < eps) {
			break;
		}

		prev_c = c;
	}

	Eigen::VectorXd f = yy - x;
	f = (f - B * solverA.solve(f)).eval();

	std::vector<double> f_result(f.size()), x_result(x.size());

	Eigen::VectorXd::Map(f_result.data(), f.size()) = f;
	Eigen::VectorXd::Map(x_result.data(), x.size()) = x;

	return std::tuple< std::vector<double>, std::vector<double> >(f_result, x_result);
}

//
// Implementation of BeadsExpandBoundaries() function
//
const std::vector<double> BeadsExpandBoundaries(const std::vector<double> & y, const unsigned int n)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BeadsExpandBoundaries(): the length of y is zero.");
	}

	std::vector<double> result(y.size() + 2 * n);

	auto it1 = result.begin();
	auto it2 = result.rbegin();
	for(unsigned int i = 0; i < n; i++, it1++, it2++) {
		double x = (double)i / (n - 1);
		double t = std::sqrt(x);
		*it1 = y.front() * t;
		*it2 = y.back() * t;
	}

	std::copy(y.begin(), y.end(), result.begin() + n);

	return result;
}

//
// Implementation of BeadsTrimBoundaries() function
//
const std::vector<double> BeadsTrimBoundaries(const std::vector<double> & y, const unsigned int n)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BeadsTrimBoundaries(): the length of y is zero.");
	}

	if(y.size() < 2 * n) {
		throw std::invalid_argument("BeadsTrimBoundaries(): the length of y is too short.");
	}

	std::vector<double> result(y.size() - 2 * n);
	auto it = y.begin() + n;

	std::copy(it, it + result.size(), result.begin());

	return result;
}

}; // namespace sablib
