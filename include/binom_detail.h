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

// binom_ln_prob_internal(), which performed basic log-factorial arithmetic,
// has been replaced by the smarter binom_ln_prob_loader() in special_func.h .


// Ok to draw this line anywhere <= 2^39 (see old binom_ln_prob_internal error
// analysis).  I've set this to 2^36 since that's roughly where I could no
// longer easily find 1 ULP deviations from the MPFR-based pbinom()
// implementation.
HEADER_INLINE uint32_t use_tdr_for_binom_lnprob(int64_t obs_tot) {
  return (obs_tot >= (1LL << 36));
}

void BinomMassMultiPPrecomp(double k, double n, dd_real* stirlerr_ddr_ptr, dd_real* half_lf_ddr_ptr);

double BinomMassJustP(double k, double n, double p, dd_real stirlerr_ddr, dd_real half_lf_ddr, uint32_t logp);

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

// n >= 2^52, min(obs_k, n - obs_k) <= 2048
double PbinomHugeTail(double obs_k, double n, td_real p_tdr, uint32_t complement, int32_t midp, uint32_t logp);

// Returns binomial distribution tail-sum when p or q is extremely small (can
// be zero).
double PbinomExtremeSuccP(double obs_k, double n, td_real p_tdr, uint32_t complement, int32_t midp, uint32_t logp);

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

dd_real binom_ltail_lik_bfrac_ddr(double obs_k, double n, dd_real p_ddr, dd_real q_ddr);

// For BinomTwoSidedP(), extends succ_odds_ratio_tdr and the incomplete one of
// {p_tdr, q_tdr} when they have only been evaluated to dd_real precision so
// far.  Does nothing if all three have already been fully evaluated.
void materialize_oddsratio_p_q_tdr(uint32_t succ_flipped, td_real* p_tdr_ptr, td_real* q_tdr_ptr, td_real* succ_odds_ratio_tdr_ptr);


#ifdef __cplusplus
}
#endif

#endif  // __BINOM_DETAIL_H__
