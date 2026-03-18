/**
 * @file beads.h
 * @brief Baseline estimation and subtraction using Baseline Estimation And Denoising using Sparsity(BEADS)
 * @author Izadori
 * @details
 * References:
 * @li Ning, X.; Selesnick, I. W.; Duval, L. Chromatogram baseline estimation and denoising using sparsity. Chemometrics and Intelligent Laboratory Systems, 2014, 139, 156-167.
 * @li [Duval's MATLAB implemetation](http://www.mathworks.com/matlabcentral/fileexchange/49974-beads--baseline-estimation-and-denoising-w--sparsity--chromatogram-signals-)
 * @li [Kotaro Saito's pybeads implementation](https://github.com/skotaro/pybeads/)
 */

#include <tuple>
#include <vector>

#ifndef __SABLIB_BEADS_H__
#define __SABLIB_BEADS_H__

namespace sablib {

/**
 * @brief Penalty types for the BEADS algorithm.(see Table 1 in Duval's paper)
 */
enum class BeadsPenalty
{
	L1_v1, /**< Penalty type 1, using phi-B(x) function and it's derivative. */
	L1_v2  /**< Penalty type 2, using phi-C(x) function and it's derivative. */
};

/**
 * @brief Performs baseline estimation and denoising using Sparsity (BEADS).
 *
 * @param y The input data.
 * @param s Order of the derivative for baseline sparsity (typically 1 or 2).
 * @param frequency Sampling frequency of the signal.
 * @param r High-pass filter parameter (cut-off frequency relative to sampling frequency).
 * @param lambda0 Sparsity parameter for the baseline.
 * @param lambda1 Sparsity parameter for the first-order derivative of the signal.
 * @param lambda2 Sparsity parameter for the second-order derivative of the signal.
 * @param loop Maximum number of iterations.
 * @param eps Convergence threshold.
 * @param penalty Penalty type (L1_v1 or L1_v2).
 * @return A tuple containing (baseline, denoised signal).
 */
const std::tuple< std::vector<double>, std::vector<double> >
BaselineBeads(
	const std::vector<double> & y, const unsigned int s, const double frequency, const double r,
	const double lambda0, const double lambda1, const double lambda2, const unsigned int loop = 30,
	const double eps = 1e-3, const BeadsPenalty penalty = BeadsPenalty::L1_v2
);

}; // namespace sablib

#endif // _SABLIB_BEADS_H__
