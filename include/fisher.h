#ifndef __FISHER_H__
#define __FISHER_H__

// Fisher's Exact Test library, copyright (C) 2013-2026 Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_base.h"

// The main ideas behind these Fisher's Exact Test and hypergeometric-function
// implementations are:
// 1. It is expensive to evaluate the (log-)probability of a contingency table
//    from scratch, but cheap to compute the likelihood-ratio between a pair of
//    adjacent contingency tables.  Both Fisher's exact test and the
//    {p,q}hyper() functions involve summing probabilities of many adjacent
//    contingency tables; we express these sums as
//      <probability of starting table> * <multiplier>
//    and calculate the multiplier by repeatedly applying the likelihood-ratio
//    formula.
// 2. Log-probabilities drop off at a faster-than-quadratic rate as one moves
//    away from the mode.  When relative errors around 2^{-53} (float64) or
//    even 2^{-24} (float32) are acceptable, it's ok to ignore contingency
//    tables with probabilities less than roughly that multiple of the starting
//    contingency table's probability; for common large cases, we can ignore
//    most tables.
// 3. For 2x3 and larger Fisher's Exact Tests, we also take advantage of
//    combinatorial identities involving sums of entire lines/planes/etc. of
//    contingency table probabilities.
// 4. With the help of the QD high-precision library, we can evaluate the
//    log-probability of a single contingency table to better-than-float64
//    precision, and accumulate partial sums with better-than-float64 precision
//    when worthwhile.  This lets us achieve much higher accuracy *and* much
//    higher speed than widely-used implementations at the same time!
// (Implementations of the Mehta/Patel network algorithm get similar mileage
// from ideas #1 and #3, but take less advantage of #2 and #4.)

#ifdef __cplusplus
namespace plink2 {
#endif

double Fisher22TwoSidedP(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, int32_t midp, uint32_t logp);

double Fisher22OddsRatio(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22);

double Fisher22OddsRatioQuantileMatch(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double target_p);

// On entry, (*lowp, *highp) should be the target quantiles.
// On exit, they are set to the odds ratios corresponding to those quantiles.
HEADER_INLINE double Fisher22OddsRatioCI(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double* lowp, double* highp) {
  // R's CI.
  *lowp = 1.0 / Fisher22OddsRatioQuantileMatch(obs_m12, obs_m11, obs_m22, obs_m21, *lowp);
  *highp = Fisher22OddsRatioQuantileMatch(obs_m11, obs_m12, obs_m21, obs_m22, 1.0 - *highp);
}

double Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp);

// Probable todos:
// - Function implementing 2xk for k>2, recursion can just follow how
//   Fisher23LnP builds on 2x2
// - Function implementing jxk for 3<=j<=k; recursive function can take a
//   parameter specifying the column number in the bottom row to iterate over
//   (bottom-row values further to the right are locked during the function
//   call).  I expect the function would rarely be worth using over simulation
//   when j>3 (complexity is roughly O(n^{(j-1)(k-1)/2})), but O(n^2) for j=k=3
//   seems likely to be a good deal.

#ifdef __cplusplus
}
#endif

#endif  // __FISHER_H__
