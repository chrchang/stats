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

#include "binom_detail.h"

#include <math.h>

#include "plink2_float.h"
#include "special_func.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Currently assumes k < n.  Should always have <1 ULP error; and
// ddr_exp(result) also has <1 ULP error when it isn't < DBL_MIN.  Error
// analysis:
//
// * For p=0.5, dd_real calculation doesn't reach absolute error > 2^{-53}
//   until n > ~2^39.
//     ddr_lfact(2^39) ~= 2^39 * log(2^39)
//                     ~= 2^39 * 39 * 0.693
//                     ~= 2^44.
//   Reviewing the operations required by ddr_lfact() (dominated by a ddr_log()
//   call, which is in turn dominated by a ddr_exp() call), I'm pretty sure
//   ddr_lfact()'s relative error < 2^{-98} for our domain (while I would be
//   surprised if it was always < 2^{-101}).  We need to add one number near
//   ddr_lfact(n) (or two numbers near ddr_lfact(n/2)) and subtract
//   ddr_lfact(n); that translates to absolute error ~2^{-53}.
//
// * For other p, if (k ln p) + ((n-k) ln q) doesn't have significantly larger
//   magnitude than the p=0.5 case, the p=0.5 error analysis applies.  If it
//   does have significantly larger magnitude, absolute error can be larger but
//   that's fine because we no longer have heavy subtractive cancellation.
dd_real binom_ln_prob_internal(int64_t k, int64_t n, dd_real p_ddr, dd_real q_ddr) {
  const uint32_t p_is_half = ddr_is(p_ddr, 0.5);
  // strictly speaking, usually don't need to initialize this at all in
  // p_is_half case
  dd_real ln_q_ddr = _ddr_log05;
  if (!p_is_half) {
    ln_q_ddr = ddr_log_2arg(q_ddr, p_ddr);
  }
  if (k == 0) {
    return ddr_muld(ln_q_ddr, n-k);
  }
  if (!use_tdr_for_binom_lnprob(n)) {
    dd_real ddrs[5];
    ddrs[0] = ddr_lfact(n);
    ddrs[1] = ddr_negate(ddr_lfact(k));
    ddrs[2] = ddr_negate(ddr_lfact(n-k));
    if (p_is_half) {
      ddrs[3] = ddr_muld(_ddr_log05, n);
    } else {
      ddrs[3] = ddr_muld(ddr_log(p_ddr), k);
      ddrs[4] = ddr_muld(ln_q_ddr, n-k);
    }
    return ddr_sort_and_add(5 - p_is_half, ddrs);
  }
  td_real tdrs[5];
  tdrs[0] = tdr_lfact(n);
  tdrs[1] = tdr_negate(tdr_lfact(k));
  tdrs[2] = tdr_negate(tdr_lfact(n-k));
  if (p_is_half) {
    tdrs[3] = tdr_muld(_tdr_log05, n);
  } else {
    // I think this consistently squeezes out enough accuracy, even for n near
    // 2^52, so no need for this function to take p_tdr?
    //   d/dp[k ln p + (n-k) ln (1-p)]
    // = k(1/p) + (n-k)(-1/(1-p))
    // = k/p - (n-k)/(1-p)
    // so if k/p ~= (n-k)/(1-p) (i.e. we're near the mode, and log-probability
    // magnitude isn't huge), a relative error of 2^{-106} in representing p
    // translates to sufficiently-low absolute error; and if we aren't near the
    // mode, the final magnitude will be large (so higher absolute error is
    // fine).
    if (ddr_ltd(p_ddr, 0.5)) {
      tdrs[3] = tdr_muld(tdr_log(tdr_make_dd(p_ddr)), k);
      tdrs[4] = tdr_muld(tdr_log1p(tdr_make_dd(ddr_negate(p_ddr))), n-k);
    } else {
      tdrs[3] = tdr_muld(tdr_log1p(tdr_make_dd(ddr_negate(q_ddr))), k);
      tdrs[4] = tdr_muld(tdr_log(tdr_make_dd(q_ddr)), n-k);
    }
  }
  return ddr_make_td(tdr_sort_and_add(5 - p_is_half, tdrs));
}

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
intptr_t BinomCompare(int64_t obs_succ, int64_t obs_tot, td_real succ_odds_ratio_tdr, int64_t succ, td_real* starting_lnprobv_tdr_ptr, td_real* ln_odds_ratio_tdr_ptr, double* dbl_ptr) {
  // Binomial probability is
  //
  //        n!        k        n-k
  //     --------- * p  * (1-p)
  //     k! (n-k)!
  //
  //                        k
  //        n!       [  p  ]         n
  //   = --------- * [ --- ]  * (1-p)
  //     k! (n-k)!   [ 1-p ]
  //
  // where k = # of successes and p is the expected success rate.
  //
  // Thus, the likelihood ratio of interest is
  //
  //   obs_succ! (obs_tot - obs_succ)!                  succ - obs_succ
  //   ------------------------------- * succ_odds_ratio
  //       succ! (obs_tot - succ)!

  uint64_t numer_factorial_args[2];
  numer_factorial_args[0] = obs_succ;
  numer_factorial_args[1] = obs_tot - obs_succ;
  uint64_t denom_factorial_args[2];
  denom_factorial_args[0] = succ;
  denom_factorial_args[1] = obs_tot - succ;
  return CompareFactorialProducts(2, succ_odds_ratio_tdr, succ - obs_succ, obs_succ, numer_factorial_args, denom_factorial_args, starting_lnprobv_tdr_ptr, ln_odds_ratio_tdr_ptr, dbl_ptr);
}

double PbinomExtremeSuccP(int64_t obs_k, int64_t n, td_real p_tdr, uint32_t complement, int32_t midp, uint32_t logp) {
  // Need to be careful about underflow, but this case is otherwise
  // straightforward since the pmf is a good-enough approximation of either the
  // cdf or ccdf.
  dd_real p_ddr = ddr_make_td(p_tdr);
  dd_real q_ddr = ddr_make_td(tdr_negate(tdr_addd(p_tdr, -1.0)));
  if (complement) {
    obs_k = n - obs_k - (!midp);
    if (obs_k < 0) {
      return logp? (0.0 / 0.0) : 0.0;
    }
    swap_ddr(&p_ddr, &q_ddr);
  }
  const double n_d = n;
  if (p_ddr.x[0] < 0.5) {
    // pmf(0) = q^n
    // pmf(1) = p * q^{n-1} * n ~= np
    // pmf(2) = p^2 * q^{n-2} * n(n-1)/2 ~= p^2 * (n(n-1)/2)
    // pmf(3 or greater) = underflow
    if (midp && (obs_k == 0)) {
      // Difference from 0.5 is too small to matter.
      return logp? -kLn2 : 0.5;
    }
    if (!logp) {
      // Difference from 1 only representable if we're in log-space.
      return 1.0;
    }
    // log(cdf(obs_k)) ~= -ccdf(obs_k) ~= -pmf(obs_k + 1)
    // log(cdf(obs_k) - 0.5 * pmf(obs_k)) ~= -0.5 * pmf(obs_k)
    const int64_t pmf_arg = obs_k + (!midp);
    if (pmf_arg > 2) {
      return 0.0;
    }
    double pmf_val;
    if (pmf_arg == 1) {
      pmf_val = ddr_muld(p_ddr, n_d).x[0];
    } else {
      // pmf_arg = 2
      // p^2 guaranteed to underflow, but p * (n(n-1)/2) * p might not.
      pmf_val = ddr_mul(ddr_muld(p_ddr, n_d * (n_d - 1) * 0.5), p_ddr).x[0];
    }
    return (0.5 * midp - 1) * pmf_val;
  }
  const int64_t nmk = n - obs_k;
  if (nmk == 0) {
    if (!midp) {
      return logp? 0 : 1;
    }
    return logp? -kLn2 : 0.5;
  }
  if (ddr_is_zero(q_ddr)) {
    return logp? (0.0 / 0.0) : 0.0;
  }
  if (nmk == 1) {
    const double pval = (1 - 0.5 * midp) * ddr_muld(q_ddr, n_d).x[0];
    return logp? log(pval) : pval;
  }
  if (!logp) {
    if (nmk > 2) {
      return 0;
    }
    // this may underflow, so don't conditionally take log
    return (1 - 0.5 * midp) * ddr_mul(ddr_muld(q_ddr, n_d * (n_d - 1) * 0.5), q_ddr).x[0];
  }
  // log(cdf(obs_k)) ~= log(pmf(obs_k))
  //                  = log(p^{obs_k} q^{nmk} (n choose nmk))
  //                 ~= log(q^{nmk} (n choose nmk))
  dd_real retval_ddr;
  if (use_tdr_for_binom_lnprob(n)) {
    dd_real ddrs[4];
    ddrs[0] = ddr_lfact(n_d);
    ddrs[1] = ddr_negate(ddr_lfact(obs_k));
    ddrs[2] = ddr_negate(ddr_lfact(nmk));
    ddrs[3] = ddr_muld(ddr_log(q_ddr), nmk);
    retval_ddr = ddr_sort_and_add(4, ddrs);
  } else {
    td_real tdrs[4];
    tdrs[0] = tdr_lfact(n_d);
    tdrs[1] = tdr_negate(tdr_lfact(obs_k));
    tdrs[2] = tdr_negate(tdr_lfact(nmk));
    tdrs[3] = tdr_muld(tdr_log(tdr_make_dd(q_ddr)), nmk);
    retval_ddr = ddr_make_td(tdr_sort_and_add(4, tdrs));
  }
  if (midp) {
    retval_ddr = ddr_sub(retval_ddr, _ddr_log2);
  }
  return retval_ddr.x[0];
}

int64_t QbinomExtremeSuccP(dd_real targetp_or_lnp_ddr, int64_t n, td_real succp_tdr, uint32_t log_target) {
  // Extremely-small-but-positive succp or failp (need to be careful about
  // underflow), 0 < targetp < 1, n > 0.
  const double n_d = n;
  if (succp_tdr.x[0] < 0.5) {
    // No, this isn't accurate when log_target is false and targetp isn't very
    // close to 1.  But we correctly return 0 in that case.
    const dd_real neg_target_lnp_ddr = ddr_negate(log_target? targetp_or_lnp_ddr : ddr_subd(targetp_or_lnp_ddr, 1));
    const dd_real succp_ddr = ddr_make_td(succp_tdr);
    const dd_real pmf1_ddr = ddr_muld(succp_ddr, n_d);
    if (ddr_leq(pmf1_ddr, neg_target_lnp_ddr)) {
      return 0;
    }
    const dd_real pmf2_ddr = ddr_mul(ddr_muld(succp_ddr, n_d * (n_d - 1) * 0.5), succp_ddr);
    return 1 + ddr_gt(pmf2_ddr, neg_target_lnp_ddr);
  }
  const dd_real failp_ddr = ddr_negate(ddr_make_td(tdr_addd(succp_tdr, -1.0)));
  if ((!log_target) || (log_target && (targetp_or_lnp_ddr.x[0] > -708.0))) {
    const dd_real targetp_ddr = log_target? ddr_exp(targetp_or_lnp_ddr) : targetp_or_lnp_ddr;
    const dd_real pmf_nm1_ddr = ddr_muld(failp_ddr, n_d);
    if (ddr_lt(pmf_nm1_ddr, targetp_ddr)) {
      return n;
    }
    const dd_real pmf_nm2_ddr = ddr_mul(ddr_muld(failp_ddr, n_d * (n_d - 1) * 0.5), failp_ddr);
    return n - 1 - ddr_geq(pmf_nm2_ddr, targetp_ddr);
  }
  // Perform search, treating ll_deriv as -log(failp).
  const uint32_t use_tdr = (n >= (1LL << 39));
  td_real logq_tdr;
  if (!use_tdr) {
    logq_tdr = tdr_make_dd(ddr_log(failp_ddr));
  } else {
    logq_tdr = tdr_log(tdr_make_dd(failp_ddr));
  }
  const dd_real target_lnp_ddr = targetp_or_lnp_ddr;
  double nmk = trunc(ddr_accurate_div(target_lnp_ddr, ddr_make_td(logq_tdr)).x[0]);
  if (nmk > n) {
    return 0;
  }
  dd_real cur_lnprob_ddr;
  dd_real diff_ddr;
  double k;
  while (1) {
    k = n_d - nmk;
    if (!use_tdr) {
      dd_real ddrs[4];
      ddrs[0] = ddr_lfact(n_d);
      ddrs[1] = ddr_negate(ddr_lfact(k));
      ddrs[2] = ddr_negate(ddr_lfact(nmk));
      ddrs[3] = ddr_muld(ddr_make_td(logq_tdr), nmk);
      cur_lnprob_ddr = ddr_sort_and_add(4, ddrs);
    } else {
      td_real tdrs[4];
      tdrs[0] = tdr_lfact(n_d);
      tdrs[1] = tdr_negate(tdr_lfact(k));
      tdrs[2] = tdr_negate(tdr_lfact(nmk));
      tdrs[3] = tdr_muld(logq_tdr, nmk);
      cur_lnprob_ddr = ddr_make_td(tdr_sort_and_add(4, tdrs));
    }
    diff_ddr = ddr_sub(cur_lnprob_ddr, target_lnp_ddr);
    if (fabs(diff_ddr.x[0]) < 16384.0) {
      break;
    }
    const double k_incr = trunc(diff_ddr.x[0] / logq_tdr.x[0]);
    nmk -= k_incr;
    if (nmk < 0) {
      nmk = 0;
    } else if (nmk > n_d) {
      nmk = n_d;
    }
  }
  const dd_real logq_ddr = ddr_make_td(logq_tdr);
  if (diff_ddr.x[0] < 0) {
    do {
      k += 1;
      diff_ddr = ddr_addd(ddr_sub(diff_ddr, logq_ddr), log(nmk / k));
      nmk -= 1;
    } while (diff_ddr.x[0] < 0.0);
    return S_CAST(int64_t, k);
  }
  do {
    nmk += 1;
    diff_ddr = ddr_addd(ddr_add(diff_ddr, logq_ddr), log(k / nmk));
    k -= 1;
  } while (diff_ddr.x[0] >= 0);
  return 1 + S_CAST(int64_t, k);
}


// Useful identities:
// 1. ibeta_power_terms_d_ln() - log p(n-k) = log-probability for succ=k
// 2. ibeta_continued_fraction_recip_d() * p(n-k) = tail-prob / succ=k
/*
dd_real binom_ln_prob_approx(int64_t k, int64_t n, dd_real p_ddr, dd_real q_ddr, double* nonlog_denom_ptr) {
  if (!((n > 512) && (MINV(k+1, n-k) >= 40))) {
    *nonlog_denom_ptr = 1;
    return binom_ln_prob_internal(k, n, p_ddr, q_ddr);
  }
  double aa = k + 1;
  double bb = n - k;
  *nonlog_denom_ptr = bb * p_ddr.x[0];
  dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), ddr_muld(p_ddr, bb));
  if (ay_minus_bx_ddr.x[0] < 0.0) {
    swap_f64(&aa, &bb);
    swap_ddr(&p_ddr, &q_ddr);
    ay_minus_bx_ddr = ddr_negate(ay_minus_bx_ddr);
  }
  return ibeta_power_terms_d_ln(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr);
}
*/

double binom_ltail_lik_simple(double succ, double fail, double succ_odds_ratio, uint32_t midp) {
  double lik = 1;
  double tail_sum = 1 - midp * 0.5;
  // Iterate outward to floating-point precision limit.
  while (1) {
    fail += 1;
    lik *= succ / (succ_odds_ratio * fail);
    succ -= 1;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      return tail_sum;
    }
  }
}

double binom_tail_lik_bfrac(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t complement, uint32_t midp) {
  double aa = obs_k + 1;
  double bb = n - obs_k;
  double xx = p_ddr.x[0];
  double yy = q_ddr.x[0];
  const double p_nmk = bb * xx;
  dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), ddr_muld(p_ddr, bb));
  uint32_t inv = !complement;
  if (ay_minus_bx_ddr.x[0] < 0.0) {
    swap_f64(&aa, &bb);
    swap_f64(&xx, &yy);
    ay_minus_bx_ddr = ddr_negate(ay_minus_bx_ddr);
    inv = !inv;
  }
  return p_nmk * ibeta_continued_fraction_recip_d(aa, bb, xx, yy, ay_minus_bx_ddr, inv, midp * (1 + complement));
}

dd_real binom_ltail_lik_simple_ddr(double k, double nmk, dd_real lik_ddr, dd_real qdp_ddr, double allowed_ulp_err) {
  if (k == 0) {
    return lik_ddr;
  }
  dd_real tailsum_ddr = lik_ddr;
  // Could use geometric-series upper bound on tailsum to raise this
  // threshold.
  const double min_incr_left = allowed_ulp_err / (k * k);
  do {
    nmk += 1;
    lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
    k -= 1;
    tailsum_ddr = ddr_add(tailsum_ddr, lik_ddr);
  } while (lik_ddr.x[0] > tailsum_ddr.x[0] * min_incr_left);
  if (k > 0) {
    const double qdp = qdp_ddr.x[0];
    double lik = lik_ddr.x[0];
    double tailsum = 0.0;
    while (1) {
      nmk += 1;
      lik *= qdp * k / nmk;
      k -= 1;
      const double preadd = tailsum;
      tailsum += lik;
      if (tailsum == preadd) {
        break;
      }
    }
    tailsum_ddr = ddr_addd(tailsum_ddr, tailsum);
  }
  return tailsum_ddr;
}

dd_real binom_ltail_lik_bfrac_ddr(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr) {
  double aa = obs_k + 1;
  double bb = n - obs_k;
  const dd_real p_nmk_ddr = ddr_muld(p_ddr, bb);
  dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), ddr_muld(p_ddr, bb));
  if (ay_minus_bx_ddr.x[0] < 0.0) {
    swap_f64(&aa, &bb);
    swap_ddr(&p_ddr, &q_ddr);
    ay_minus_bx_ddr = ddr_negate(ay_minus_bx_ddr);
  }
  return ddr_accurate_div(p_nmk_ddr, ibeta_continued_fraction_ddr(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr));
}

void materialize_oddsratio_p_q_tdr(uint32_t succ_flipped, td_real* p_tdr_ptr, td_real* q_tdr_ptr, td_real* succ_odds_ratio_tdr_ptr) {
  // Currently safe to assume that either succ_odds_ratio_tdr is fully
  // evaluated, or both succ_odds_ratio_tdr and one of {p_tdr, q_tdr} needs to
  // be.
  if (succ_odds_ratio_tdr_ptr->x[2] == DBL_MAX) {
    if (!succ_flipped) {
      *q_tdr_ptr = tdr_addd(tdr_negate(*p_tdr_ptr), 1.0);
    } else {
      *p_tdr_ptr = tdr_addd(tdr_negate(*q_tdr_ptr), 1.0);
    }
    *succ_odds_ratio_tdr_ptr = tdr_accurate_div(*p_tdr_ptr, *q_tdr_ptr);
  }
}

#ifdef __cplusplus
}
#endif
