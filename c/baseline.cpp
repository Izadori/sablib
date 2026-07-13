/**
 * @file baseline.cpp
 * @brief The C interface of sablib baseline estimation functions.(implementation)
 * @author Izadori
 */

#include "baseline.h"

#include "../baseline/airpls.h"
#include "../baseline/arpls.h"
#include "../baseline/asls.h"
#include "../baseline/aspls.h"
#include "../baseline/backcor.h"
#include "../baseline/beads.h"
#include "../baseline/goldindec.h"
#include "../baseline/imodpoly.h"
#include "../baseline/modpoly.h"
#include "../baseline/polynomial.h"
#include "../baseline/psalsa.h"
#include "../baseline/sma.h"
#include "../baseline/snip.h"
#include "../baseline/spline.h"

//
// Implementation of Sablib_BaselineLinear() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineLinear(
    const SABLIB_DATA_PTR y, const unsigned int index1, const unsigned int index2
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineLinear(yy, index1, index2);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselinePolynomial() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselinePolynomial(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const unsigned int * indices_ptr, const size_t ptr_size
)
{
    std::vector<double> yy;
    std::vector<unsigned int> indices;
    yy.assign(y->data, y->data + y->size);
    indices.assign(indices_ptr, indices_ptr + ptr_size);

    auto result = sablib::BaselinePolynomial(yy, polyorder, indices);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineSpline() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineSpline(
    const SABLIB_DATA_PTR y, const unsigned int * indices_ptr, const size_t ptr_size
)
{
    std::vector<double> yy;
    std::vector<unsigned int> indices;
    yy.assign(y->data, y->data + y->size);
    indices.assign(indices_ptr, indices_ptr + ptr_size);

    auto result = sablib::BaselineSpline(yy, indices);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineSMA() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineSMA(
    const SABLIB_DATA_PTR y, const unsigned int n, const unsigned int loop
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineSMA(yy, n, loop);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineSnip() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineSnip(
	const SABLIB_DATA_PTR y, const unsigned int m, const bool decreasing,
	const enum Sablib_SnipPreprocess preprocess, const unsigned int loop
)
{
    sablib::SnipPreprocess pp;

    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    switch(preprocess) {
    case LL:
        pp = sablib::SnipPreprocess::LL;
        break;
    case LLS:
        pp = sablib::SnipPreprocess::LLS;
        break;
    default:
        pp = sablib::SnipPreprocess::None;
        break;
    };

    auto result = sablib::BaselineSnip(yy, m, decreasing, pp, loop);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineModPoly() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineModPoly(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineModPoly(yy, polyorder, loop, eps);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineIModPoly() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineIModPoly(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const double k,
	const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineIModPoly(yy, polyorder, k, loop, eps);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineBackcor() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineBackcor(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const enum Sablib_BackcorFunc func,
	const double s, const double alpha, const unsigned int loop, const double eps
)
{
    sablib::BackcorFunc f;

    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    switch(func) {
    case Huber:
        f = sablib::BackcorFunc::Huber;
        break;
    case AHuber:
        f = sablib::BackcorFunc::AHuber;
        break;
    case TQuad:
        f = sablib::BackcorFunc::TQuad;
        break;
    case Indec:
        f = sablib::BackcorFunc::Indec;
        break;
    case AIndec:
        f = sablib::BackcorFunc::AIndec;
        break;
    default:
        f = sablib::BackcorFunc::ATQuad;
        break;
    }

    auto result = sablib::BaselineBackcor(yy, polyorder, f, s, alpha, loop, eps);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineGoldindec() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineGoldindec(
	const SABLIB_DATA_PTR y, const unsigned int polyorder, const double peak_ratio,
	const double alpha, const unsigned int loop, const double eps,
	const unsigned int loop_legend, const double eps_legend, const double eps_s
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineGoldindec(
        yy, polyorder, peak_ratio, alpha, loop, eps, loop_legend, eps_legend, eps_s
    );

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineAsLS() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAsLS(
	const SABLIB_DATA_PTR y, const double lambda, const double p, const unsigned int s,
	const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineAsLS(yy, lambda, p, s, loop, eps);

    SABLIB_BASELINE_DATA_PTR ptr = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), ptr->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), ptr->corrected.data);

    return ptr;
}

//
// Implementation of Sablib_BaselineAirPLS() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAirPLS(
	const SABLIB_DATA_PTR y, const double lambda, const unsigned int s,
	const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineAirPLS(yy, lambda, s, loop, eps);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineArPLS() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineArPLS(
	const SABLIB_DATA_PTR y, const double lambda, const unsigned int s,
	const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineArPLS(yy, lambda, s, loop, eps);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselineAsPLS() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineAsPLS(
	const SABLIB_DATA_PTR y, const double lambda, const double k,
    const unsigned int s, const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselineAsPLS(yy, lambda, k, s, loop, eps);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BaselinePsalsa() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselinePsalsa(
	const SABLIB_DATA_PTR y, const double lambda, const double p, const double k,
	const unsigned int s, const unsigned int loop, const double eps
)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BaselinePsalsa(yy, lambda, p, k, s, loop, eps);

    SABLIB_BASELINE_DATA_PTR ptr = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), ptr->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), ptr->corrected.data);

    return ptr;
}

//
// Implementation of Sablib_BaselineBeads() function.
//
const SABLIB_BASELINE_DATA_PTR Sablib_BaselineBeads(
	const SABLIB_DATA_PTR y, const unsigned int s, const double frequency, const double r,
	const double lambda0, const double lambda1, const double lambda2, const unsigned int loop,
	const double eps, const enum Sablib_BeadsPenalty penalty
)
{
    sablib::BeadsPenalty pen;

    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    switch(penalty) {
    case L1_v1:
        pen = sablib::BeadsPenalty::L1_v1;
        break;
    default:
        pen = sablib::BeadsPenalty::L1_v2;
        break;
    }

    auto result = sablib::BaselineBeads(yy, s, frequency, r, lambda0, lambda1, lambda2, loop, eps, pen);

    SABLIB_BASELINE_DATA_PTR p = AllocSablibBaselineData(result.baseline.size());
    std::copy(result.baseline.begin(), result.baseline.end(), p->baseline.data);
    std::copy(result.corrected.begin(), result.corrected.end(), p->corrected.data);

    return p;
}

//
// Implementation of Sablib_BeadsExpandBoundaries() function.
//
const SABLIB_DATA_PTR Sablib_BeadsExpandBoundaries(const SABLIB_DATA_PTR y, const unsigned int n)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BeadsExpandBoundaries(yy, n);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}

//
// Implementation of Sablib_BeadsTrimBoundaries() function.
//
const SABLIB_DATA_PTR Sablib_BeadsTrimBoundaries(const SABLIB_DATA_PTR y, const unsigned int n)
{
    std::vector<double> yy;
    yy.assign(y->data, y->data + y->size);

    auto result = sablib::BeadsTrimBoundaries(yy, n);

    SABLIB_DATA_PTR p = AllocSablibData(result.size());
    std::copy(result.begin(), result.end(), p->data);

    return p;
}
