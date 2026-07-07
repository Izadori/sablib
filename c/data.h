/**
 * @file data.h
 * @brief Definition of data structures for the C interface of sablib, creation and release of data.
 * @author Izadori
 */

#ifndef __SABLIB_C_DATA_H__
#define __SABLIB_C_DATA_H__

#ifdef __cplusplus
#include <cstddef>  // for size_t
#else
#include <stddef.h> // for size_t
#endif // __cplusplus

/**
 * @brief Definition of sablib data structure.
 */
typedef struct _stSablibData
{
    size_t size;   /**< Data size. */
    double * data; /**< Data body. */
} SABLIB_DATA, *SABLIB_DATA_PTR;

/**
 * @brief Definition of the data structure for the return value of sablib baseline functions.
 */
typedef struct _stSablibBaselineData
{
    SABLIB_DATA baseline;  /**< Estimated baseline. */
    SABLIB_DATA corrected; /**< Corrected data (data with baseline subtracted). */
} SABLIB_BASELINE_DATA, *SABLIB_BASELINE_DATA_PTR;

#include "sablib_export.h"

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

/**
 * @brief Allocates a SABLIB_DATA structure.
 * @param size Size of the data to allocate.
 * @return A pointer to the allocated SABLIB_DATA structure, or NULL if allocation fails.
 */
SABLIB_EXPORT const SABLIB_DATA_PTR AllocSablibData(const size_t size);

/**
 * @brief Frees a SABLIB_DATA structure and its data body.
 * @param data Pointer to the SABLIB_DATA structure to free.
 */
SABLIB_EXPORT void FreeSablibData(SABLIB_DATA_PTR data);

/**
 * @brief Allocates a SABLIB_BASELINE_DATA structure.
 * @param size Size of the data to allocate for baseline and corrected data.
 * @return A pointer to the allocated SABLIB_BASELINE_DATA structure, or NULL if allocation fails.
 */
SABLIB_EXPORT const SABLIB_BASELINE_DATA_PTR AllocSablibBaselineData(const size_t size);

/**
 * @brief Frees a SABLIB_BASELINE_DATA structure and its members.
 * @param data Pointer to the SABLIB_BASELINE_DATA structure to free.
 */
SABLIB_EXPORT void FreeSablibBaselineData(SABLIB_BASELINE_DATA_PTR data);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // __SABLIB_C_DATA_H__
