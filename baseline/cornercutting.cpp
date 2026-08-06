/**
 * @file cornercutting.cpp
 * @brief Baseline estimation using Corner-Cutting method(implementation)
 * @author Izadori
 */

#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>

#include "cornercutting.h"

namespace sablib {

namespace {

const std::vector<double> BezierInterpolate(
	const std::vector<double> & y, const std::vector<unsigned int> & index
)
{
	struct point {
		double x;
		double y;
	};

	std::vector<point> p;
	std::vector<double> z(y.size());

	p.push_back({(double)index[0], y[index[0]]});

	for(unsigned int i = 1; i < index.size() - 2; i++) {
		p.push_back({(double)index[i], y[index[i]]});
		p.push_back({(double)(index[i] + index[i + 1]) / 2, (y[index[i]] + y[index[i + 1]]) / 2});
	}

	p.push_back({(double)index[index.size() - 2], y[index[index.size() - 2]]});
	p.push_back({(double)index[index.size() - 1], y[index[index.size() - 1]]});

	for(unsigned int i = 0, j = 0; j < p.size() - 2; j += 2) {
		double a = p[j].x - 2 * p[j + 1].x + p[j + 2].x;
		double b = -2 * p[j].x + 2 * p[j + 1].x;

		while((double)i <= p[j + 2].x) {
			double c = p[j].x - (double)i;
			double t = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);

			z[i] = (1 - t) * (1 - t) * p[j].y + 2 * (1 - t) * t * p[j + 1].y + t * t * p[j + 2].y;
			i++;
		}
	}

	return z;
}

}; // unnamed namespace

//
// Implementation of BaselineCornerCutting() function
//
const BaselineResult
BaselineCornerCutting(
	const std::vector<double> & y, const unsigned int loop
)
{
	if(y.size() == 0) {
		throw std::invalid_argument("BaselineCornerCutting(): the length of y is zero.");
	}

	if(loop == 0) {
		throw std::invalid_argument("BaselineCornerCutting(): loop is zero.");
	}

	std::vector<unsigned int> t(y.size(), 0);
	std::vector<double> er_list;
	double a_old = 0;
	unsigned int marked_size_old = y.size();

	for(unsigned int j = 0; j < y.size() - 1; j++) {
		a_old += (y[j + 1] + y[j]);
	}

	a_old /= 2;

	for(unsigned int i = 0; i < loop; i++) {
		std::vector<unsigned int> p;

		for(unsigned int j = 0; j < t.size(); j++) {
			if(t[j] == i) {
				p.push_back(j);
			}
		}

		std::vector<unsigned int> marked;
		unsigned int marked_size_new;
		unsigned int elmininated_count = 0;

		marked.push_back(p.front());
		t[p.front()]++;

		for(unsigned int j = 1; j < p.size() - 1; j++) {
			double c = (y[p[j + 1]] - y[p[j - 1]]) / (p[j + 1] - p[j - 1]) * (p[j] - p[j - 1]) + y[p[j - 1]];

			if(y[p[j]] < c) {
				marked.push_back(p[j]);
				t[p[j]]++;
			}
		}

		marked.push_back(p.back());
		t[p.back()]++;

		marked_size_new = marked.size();
		elmininated_count = marked_size_old - marked_size_new;

		if(marked_size_new == 0 || elmininated_count == 0) {
			break;
		}

		// Note:
		// A straightforward reading of the Terminal Condition in Liu's paper suggests
		// that the area should be calculated using the remaining points. However,
		// Algorithm 2 can also be interpreted as calculating the area using the
		// removed points. Furthermore, the definition of ER is also different.
		// In this implementation, we give priority to the description in the
		// Terminal Condition and calculate the area using the remaining points.
		double a_new = 0;

		for(unsigned int j = 0; j < marked.size() - 1; j++) {
			a_new += ((double)marked[j + 1] - (double)marked[j]) * (y[marked[j + 1]] + y[marked[j]]);
		}

		a_new /= 2;

		// Calculate ER value.
		double er = (a_old - a_new) / elmininated_count;
		er_list.push_back(er);
		a_old = a_new;
		marked_size_old = marked_size_new;
	}

	// Search for maximum ER value
	unsigned int max_er_index = std::max_element(er_list.begin(), er_list.end()) - er_list.begin();

	// Make index list
	std::vector<unsigned int> index;

	// Note:
	// We search for the maximum ER value, but note that the points used to
	// calculate the ER value (i.e., the points used to construct the Bézier
	// curve) are those before T is incremented.
	// Therefore, the point corresponding to the maximum ER value can only be
	// identified by using the ER value from the previous iteration.
	for(unsigned int i = 0; i < t.size(); i++) {
		if(t[i] >= max_er_index - 1) {
			index.push_back(i);
		}
	}

	// Corner-Cutting with Bezier technique
	std::vector<double> z = BezierInterpolate(y, index);

	BaselineResult result;
	result.baseline.resize(z.size());
	result.corrected.resize(z.size());

	for(unsigned int i = 0; i < z.size(); i++) {
		result.baseline[i] = z[i];
		result.corrected[i] = y[i] - z[i];
	}

	return result;
}

}; // namespace sablib
