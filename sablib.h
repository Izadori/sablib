/**
 * @file sablib.h
 * @brief Header file of sablib, a C++ library for signal smoothing and baseline estimation
 * @author Izadori
 */

#ifndef __SABLIB_H__
#define __SABLIB_H__

#include "smoothing/moving_average.h"
#include "smoothing/savitzky_golay.h"
#include "smoothing/whittaker.h"
#include "smoothing/pspline.h"

#include "baseline/sma.h"
#include "baseline/polynomial.h"
#include "baseline/spline.h"
#include "baseline/modpoly.h"
#include "baseline/asls.h"
#include "baseline/airpls.h"
#include "baseline/arpls.h"
#include "baseline/beads.h"

#endif // __SABLIB_H__
