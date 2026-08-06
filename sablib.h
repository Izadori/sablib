/**
 * @file sablib.h
 * @brief Header file of sablib, a C++ library for signal smoothing and baseline estimation
 * @author Izadori
 */

#ifndef __SABLIB_H__
#define __SABLIB_H__

#include "smoothing/moving_average.h"
#include "smoothing/moving_median.h"
#include "smoothing/savitzky_golay.h"
#include "smoothing/whittaker.h"
#include "smoothing/pspline.h"

#include "baseline/sma.h"
#include "baseline/snip.h"
#include "baseline/polynomial.h"
#include "baseline/spline.h"
#include "baseline/modpoly.h"
#include "baseline/imodpoly.h"
#include "baseline/backcor.h"
#include "baseline/goldindec.h"
#include "baseline/asls.h"
#include "baseline/iasls.h"
#include "baseline/airpls.h"
#include "baseline/aispls.h"
#include "baseline/arpls.h"
#include "baseline/aspls.h"
#include "baseline/psalsa.h"
#include "baseline/beads.h"
#include "baseline/rubberband.h"
#include "baseline/kajfoszkwiatek.h"
#include "baseline/cornercutting.h"

#endif // __SABLIB_H__
