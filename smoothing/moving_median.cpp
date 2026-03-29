/**
 * @file moving_median.cpp
 * @brief Smoothing using moving median(implementation)
 * @author Izadori
 */

#include <algorithm>
#include <set>

#include "moving_median.h"

namespace sablib {

namespace {

//
// Implementation of SlidingMedian class
//
class SlidingMedian
{
public:
	SlidingMedian(const unsigned int n) : window_size(2 * n + 1) {};

	void Add(const double value)
	{
		if(low.empty() || value <= *low.rbegin()) {
			low.insert(value);
		}
		else {
			high.insert(value);
		}

		Rebalance();
	}

	void Remove(const double value)
	{
		auto it = low.find(value);

		if(it != low.end()) {
			low.erase(it);
		}
		else {
			it = high.find(value);

			if(it != high.end()) {
				high.erase(it);
			}
		}

		Rebalance();
	}

	double Median() const
	{
		if(window_size % 2 == 1) {
			return *low.rbegin();
		}
		else {
			return (*low.rbegin() + *high.begin()) / 2;
		}
	}

private:
	unsigned int window_size;
	// low : lower value group
	// high: higher value group
	std::multiset<double> high, low;

	void Rebalance()
	{
		if(low.size() > high.size() + 1) {
			high.insert(*low.rbegin());
			low.erase(low.find(*low.rbegin()));
		}
		else if(high.size() > low.size()) {
			low.insert(*high.begin());
			high.erase(high.begin());
		}
	}
};

}; // unnamed namespace

//
// Implementation of MovingMedian() function
//
const std::vector<double> MovingMedian(const std::vector<double> & y, const unsigned int n)
{
	if(y.size() == 0) {
		throw std::invalid_argument("MovingMedian(): the length of y is zero.");
	}

	if(n == 0) {
		throw std::invalid_argument("MovingMedian(): n is zero.");
	}

	unsigned int points = 2 * n + 1;
	SlidingMedian smed(n);
	std::vector<double> yy(y.size() + 2 * n), result;

	std::fill_n(yy.begin(), n, y.front());
	std::copy(y.begin(), y.end(), yy.begin() + n);
	std::fill_n(yy.begin() + n + y.size(), n, y.back());

	result.reserve(y.size());

	for(unsigned int i = 0; i < yy.size(); i++) {
		smed.Add(yy[i]);

		if(i >= points) {
			smed.Remove(yy[i - points]);
		}

		if(i > points) {
			result.emplace_back(smed.Median());
		}
	}

	return result;
}

}; // namespace sablib
