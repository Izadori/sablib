/**
 * @file data.cpp
 * @brief creation and release of data for the C interface of sablib.(implementation)
 * @author Izadori
 */

#include "data.h"

//
// Implementation of AllocSablibData() function
//
const SABLIB_DATA_PTR AllocSablibData(const size_t size)
{
    SABLIB_DATA_PTR p = new SABLIB_DATA;
    p->size = size;
    p->data = new double[size];

    return p;
}

//
// Implementation of FreeSablibData() function
//
void FreeSablibData(SABLIB_DATA_PTR data)
{
    if(data != nullptr) {
        delete data->data;
        delete data;
    }

    data = nullptr;
}

//
// Implementation of AllocSablibBaselineData() function
//
const SABLIB_BASELINE_DATA_PTR AllocSablibBaselineData(const size_t size)
{
    SABLIB_BASELINE_DATA_PTR p = new SABLIB_BASELINE_DATA;
    p->baseline.size = size;
    p->baseline.data = new double[size];
    p->corrected.size = size;
    p->corrected.data = new double[size];

    return p;
}

//
// Implementation of FreeSablibBaselineData() function
//
void FreeSablibBaselineData(SABLIB_BASELINE_DATA_PTR data)
{
    if(data != nullptr) {
        delete data->baseline.data;
        delete data->corrected.data;
        delete data;
    }

    data = nullptr;
}
