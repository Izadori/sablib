/**
 * @file smoothing.cpp
 * @brief The C interface of sablib smoothing functions.(implementation)
 * @author Izadori
 */

#include <algorithm>
#include <vector>

#include "smoothing.h"

#include "../smoothing/moving_average.h"
#include "../smoothing/moving_median.h"
#include "../smoothing/pspline.h"
#include "../smoothing/savitzky_golay.h"
#include "../smoothing/whittaker.h"

//
// Implementation of Sablib_MovingAverage() function.
//
const SABLIB_DATA_PTR Sablib_MovingAverage(const SABLIB_DATA_PTR y, const unsigned int n)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::MovingAverage(yy, n);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_WeightedMovingAverage() function.
//
const SABLIB_DATA_PTR Sablib_WeightedMovingAverage(const SABLIB_DATA_PTR y, const SABLIB_DATA_PTR w)
{
    std::vector<double> yy, ww;
    yy.assign(y->data, y->data + y->size);
    ww.assign(w->data, w->data + w->size);

    auto result = sablib::WeightedMovingAverage(yy, ww);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_GaussianKernel() function.
//
const SABLIB_DATA_PTR Sablib_GaussianKernel(const unsigned int n, const double sigma)
{
    auto result = sablib::GaussianKernel(n, sigma);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_GaussianFilter() function.
//
const SABLIB_DATA_PTR Sablib_GaussianFilter(const SABLIB_DATA_PTR y, const unsigned int n, const double sigma)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::GaussianFilter(yy, n, sigma);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_MovingMedian() function.
//
const SABLIB_DATA_PTR Sablib_MovingMedian(const SABLIB_DATA_PTR y, const unsigned int n)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::MovingMedian(yy, n);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_PSpline() function.
//
const SABLIB_DATA_PTR Sablib_PSpline(
	const SABLIB_DATA_PTR y, const unsigned int knots_num,
	const unsigned int degree, const unsigned int s, const double lambda
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::PSpline(yy, knots_num, degree, s, lambda);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_SavitzkyGolayCoefficients() function.
//
const SABLIB_DATA_PTR Sablib_SavitzkyGolayCoefficients(
	const unsigned int n, const unsigned int polyorder, const unsigned derive, const double delta
)
{
    auto result = sablib::SavitzkyGolayCoefficients(n, polyorder, derive, delta);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_SavitzkyGolay() function.
//
const SABLIB_DATA_PTR Sablib_SavitzkyGolay(
	const SABLIB_DATA_PTR y,
	const unsigned int n, const unsigned int polyorder, const unsigned derive, const double delta
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::SavitzkyGolay(yy, n, polyorder, derive, delta);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_WeightedWhittaker() function.
//
const SABLIB_DATA_PTR Sablib_WeightedWhittaker(
	const SABLIB_DATA_PTR y, const SABLIB_DATA_PTR w,
	const double lambda, const unsigned int s
)
{
    std::vector<double> yy, ww;
    yy.assign(y->data, y->data + y->size);
    ww.assign(w->data, w->data + w->size);

    auto result = sablib::Whittaker(yy, ww, lambda, s);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_Whittaker() function.
//
const SABLIB_DATA_PTR Sablib_Whittaker(
	const SABLIB_DATA_PTR y, const double lambda, const unsigned int s
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::Whittaker(yy, lambda, s);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}
