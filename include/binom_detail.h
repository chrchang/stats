#ifndef __BINOM_DETAIL_H__
#define __BINOM_DETAIL_H__

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

// Support routines specific to the binomial distribution.

// Ok to draw this line anywhere <= 2^39 (see binom_ln_prob_internal error
// analysis).  I've set this to 2^36 since that's roughly where I could no
// longer easily find 1 ULP deviations from the MPFR-based pbinom
// implementation.
HEADER_INLINE uint32_t use_tdr_for_binom_lnprob(int64_t obs_tot) {
  return (obs_tot >= (1LL << 36));
}

// Should always have <1 ULP error; and ddr_exp(result) also has <1 ULP error
// when it isn't < DBL_MIN.
//
// This implementation is a bit slow, but it's relatively simple and reliable.
// (ibeta_power_terms_d_ln() trades off a tiny bit of accuracy for a
// significant speed improvement.)
dd_real binom_ln_prob_internal(int64_t k, int64_t n, dd_real p_ddr, dd_real q_ddr);

void BinomMassMultiKPrecomp(int64_t n, td_real p_tdr, uint32_t* p_is_half_ptr, td_real* lfact_n_tdr_ptr, td_real* lnp_tdr_ptr, td_real* lnq_tdr_ptr);

double BinomMassJustK(int64_t k, int64_t n, uint32_t p_is_half, const td_real lfact_n_tdr, const td_real lnp_tdr, const td_real lnq_tdr, uint32_t logp);

void BinomMassMultiPPrecomp(int64_t k, int64_t n, td_real* lfact_n_tdr_ptr, td_real* neg_lfact_k_tdr_ptr, td_real* neg_lfact_nmk_tdr_ptr);

double BinomMassJustP(td_real p_tdr, int64_t k, int64_t n, const td_real lfact_n_tdr, const td_real neg_lfact_k_tdr, const td_real neg_lfact_nmk_tdr, uint32_t logp);

// - succ_odds_ratio_tdr must be p/(1-p), where p is the expected success rate.
//
// - starting_lnprobv_tdr is expected to either be initialized to
//     log(succ_odds_ratio^obs_succ / (obs_succ! (obs_tot - obs_succ)!)),
//   possibly with x[2] initialized to DBL_MAX to indicate that the calculation
//   has only been executed to dd_real precision so far; or have x[0]
//   initialized to DBL_MAX to indicate that the calculation hasn't yet
//   happened at all.  On return, the value may be refined.
//
// - ln_odds_ratio_tdr is expected to either be initialized to log(odds_ratio),
//   or have x[0] initialized to DBL_MAX, etc.
//
// - Return value is positive if succ has higher probability than obs_succ, 0
//   if identical probability, and negative if lower probability.
intptr_t BinomCompare(int64_t obs_succ, int64_t obs_tot, td_real succ_odds_ratio_tdr, int64_t succ, td_real* starting_lnprobv_tdr_ptr, td_real* ln_odds_ratio_tdr_ptr, double* dbl_ptr);

// Returns binomial distribution tail-sum when p or q is extremely small (can
// be zero).
double PbinomExtremeSuccP(int64_t obs_k, int64_t n, td_real p_tdr, uint32_t complement, int32_t midp, uint32_t logp);

// Returns smallest k for which cdf(k) >= targetp, when succp or failp is
// positive but extremely small.  Also assumes n > 0.
int64_t QbinomExtremeSuccP(dd_real targetp_or_lnp_ddr, int64_t n, td_real succp_tdr, uint32_t log_target);

// Returns binomial left-tail relative-likelihood, evaluated to ordinary
// accuracy.
//
//   pmf(succ) + pmf(succ-1) + ... + pmf(0)
//   --------------------------------------
//                 pmf(succ)
//
// with pmf(succ)/2 subtracted from the numerator when midp is true.
double binom_ltail_lik_simple(double succ, double fail, double succ_odds_ratio, uint32_t midp);

// ibeta_continued_fraction_recip_d() wrapper mirroring
// binom_ltail_lik_simple().
double binom_tail_lik_bfrac(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t complement, uint32_t midp);

// High-accuracy versions of the above.
dd_real binom_ltail_lik_simple_ddr(double k, double nmk, dd_real lik_ddr, dd_real qdp_ddr, double allowed_ulp_err);

dd_real binom_ltail_lik_bfrac_ddr(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr);

// For BinomTwoSidedP(), extends succ_odds_ratio_tdr and the incomplete one of
// {p_tdr, q_tdr} when they have only been evaluated to dd_real precision so
// far.  Does nothing if all three have already been fully evaluated.
void materialize_oddsratio_p_q_tdr(uint32_t succ_flipped, td_real* p_tdr_ptr, td_real* q_tdr_ptr, td_real* succ_odds_ratio_tdr_ptr);


#ifdef __cplusplus
}
#endif

#endif  // __BINOM_DETAIL_H__
