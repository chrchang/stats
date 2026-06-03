#ifndef __BINOM_H__
#define __BINOM_H__

// Binomial Exact Test library, copyright (C) 2013-2026 Christopher Chang.
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
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// The main ideas behind the binomial distribution and test implementations
// are:
// 1. It is expensive to evaluate the (log-)probability of an outcome from
//    scratch, but cheap to compute the likelihood-ratio between a pair of
//    adjacent outcomes.  Both the binomial test and the {p,q}binom() functions
//    involve summing probabilities of many adjacent contingency tables; we can
//    express these sums as
//      <probability of starting table> * <multiplier>
//    and calculate the multiplier by repeatedly applying the likelihood-ratio
//    formula.
// 2. Log-probabilities drop off at a faster-than-quadratic rate as one moves
//    away from the mode.  When relative errors around 2^{-53} (float64) or
//    even 2^{-24} (float32) are acceptable, it's ok to ignore outcomes with
//    probabilities less than ~that multiple of the starting outcome's
//    probability; for common large cases, we can ignore most outcomes.
// 3. (TODO: summary of DiDonato and Morris's BFRAC)
// 4. With the help of the QD high-precision library, we can evaluate the
//    log-probability of a single outcome to better-than-float64 precision, and
//    accumulate partial sums with better-than-float64 precision when
//    worthwhile.

double LnBinomCoeff(int64_t n, int64_t k);

double BinomMass(int64_t k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t logp);

double PbinomApprox(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t complement, int32_t midp, uint32_t logp);

HEADER_INLINE double BinomOneSidedP(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t succ_is_greater_alt, int32_t midp, uint32_t logp) {
  const int64_t k_decr = succ_is_greater_alt && (!midp);
  if (k_decr && (obs_k == 0)) {
    return logp? 0.0 : 1.0;
  }
  return PbinomApprox(obs_k - k_decr, n, p_ddr, q_ddr, succ_is_greater_alt, midp, logp);
}

double Pbinom(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t complement, uint32_t logp);

int64_t Qbinom(dd_real targetp_or_lnp_ddr, int64_t n, dd_real succp_ddr, uint32_t log_target);

HEADER_INLINE int64_t QbinomHalfUlp(dd_real targetp_or_lnp_ddr, int64_t n, dd_real succp_ddr, uint32_t log_target) {
  if (!ddr_is_zero(targetp_or_lnp_ddr)) {
    const double half_ulp = 0.5 * fabs(targetp_or_lnp_ddr.x[0] - prev_float64(targetp_or_lnp_ddr.x[0]));
    targetp_or_lnp_ddr = ddr_subd(targetp_or_lnp_ddr, half_ulp);
  }
  return Qbinom(targetp_or_lnp_ddr, n, succp_ddr, log_target);
}

double BinomTwoSidedP(int64_t obs_succ, int64_t obs_tot, td_real p_tdr, int32_t midp, uint32_t logp);

#ifdef __cplusplus
}
#endif

#endif  // __BINOM_H__
