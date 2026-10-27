#ifndef __FISHER_H__
#define __FISHER_H__

// Fisher's Exact Test library, copyright (C) 2013-2026 Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_base.h"
#include "plink2_highprec.h"

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

double Fisher22TwoSidedPEx(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, td_real odds_tdr, int32_t midp, uint32_t logp);

double Fisher22OddsRatio(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22);

double Fisher22OddsRatioQuantileMatch(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double target_p);

// (low_p, high_p) should be the target quantiles.
// On exit, (*low_resultp, *high_resultp) is set to the odds ratios
// corresponding to those quantiles.
HEADER_INLINE void Fisher22OddsRatioCI(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double low, double high, double* low_resultp, double* high_resultp) {
  // R's CI.
  *low_resultp = 1.0 / Fisher22OddsRatioQuantileMatch(obs_m12, obs_m11, obs_m22, obs_m21, low);
  *high_resultp = Fisher22OddsRatioQuantileMatch(obs_m11, obs_m12, obs_m21, obs_m22, 1.0 - high);
}

double Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp);

// obs[] is in row-major order: obs[0..(col_ct-1)] is the first row, and
// obs[col_ct..(2*col_ct-1)] is the second row.  Entries must be nonnegative
// and sum to less than 2^31.
// Returns log-[mid]p-value, or NaN on allocation failure.
// If there are more than 3 nonempty columns, and a cheap upper bound on the
// log-p-value is less than ln_p_cutoff, that upper bound is returned instead.
// (The bound is ln(observed table probability) + ln(bound on number of tables
// with the same margins), and since the p-value is at least the observed
// table probability, it overshoots by less than the latter term, which is
// roughly (k-1) * ln(n).  So this mainly helps when the p-value is tiny.)
// Pass -DBL_MAX to always compute the exact value.
// Complexity is roughly O(n^{(k-1)/2}) where k is the number of nonempty
// columns, so this is mainly useful for small k.
double Fisher2xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff);

// obs[] is a 3x3 table in row-major order.  Entries must be nonnegative and
// sum to less than 2^31.  Other semantics are the same as for Fisher2xkLnP().
// Complexity is roughly O(n^2).
double Fisher33LnP(const int32_t* obs, uint32_t midp, double ln_p_cutoff);

// Probable todos:
// - Larger jxk tables.  Fisher33LnP()'s implementation generalizes to 3xk,
//   but it's usually slower than the network algorithm for k > 3; merging
//   subproblems with the same remaining row totals (as the network algorithm
//   does) would probably be necessary for that to be worthwhile.

#ifdef __cplusplus
}
#endif

#endif  // __FISHER_H__
