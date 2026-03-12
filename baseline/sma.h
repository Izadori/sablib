//
// sma.h
//
// Copyright (c) 2026 Izadori
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef __SABLIB_SMA_H__
#define __SABLIB_SMA_H__

#include "../smoothing/moving_average.h"

namespace sablib {

std::vector<double> BaselineSMA(std::vector<double> & y, const unsigned int n, const unsigned int loop = 50);

}; // namespace sablib

#endif // __SABLIB_SMA_H__
