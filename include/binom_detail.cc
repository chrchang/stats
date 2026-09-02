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

// Low-level interface for multiple-p vectorized dbinom().
// (Loader's algorithm can't be sped up much in the multiple-k case.)
void BinomMassMultiPPrecomp(double k, double n, dd_real* stirlerr_ddr_ptr, dd_real* half_lf_ddr_ptr) {
  if ((k > 0) && (k < n)) {
    binom_ln_prob_loader_part1(ddr_maked(k), ddr_add2d(n, -k), ddr_maked(n), stirlerr_ddr_ptr, half_lf_ddr_ptr);
  }
}

double BinomMassJustP(double k, double n, double p, dd_real stirlerr_ddr, dd_real half_lf_ddr, uint32_t logp) {
  const dd_real n_ddr = ddr_maked(n);
  const dd_real p_ddr = ddr_maked(p);
  const dd_real q_ddr = ddr_add2d(1.0, -p);
  dd_real ln_prob_ddr;
  if (k == 0) {
    ln_prob_ddr = ddr_mul(ddr_log_extdomain_maybehalf(q_ddr), n_ddr);
  } else if (k == n) {
    ln_prob_ddr = ddr_mul(ddr_log_extdomain_maybehalf(p_ddr), n_ddr);
  } else {
    ln_prob_ddr = binom_ln_prob_loader_part2(ddr_maked(k), ddr_add2d(n, -k), n_ddr, p_ddr, q_ddr, stirlerr_ddr, half_lf_ddr);
  }
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  return ddr_exp(ln_prob_ddr).x[0];
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
//
// (Possible to do better with td_real implementation of Loader's algorithm.)
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

double PbinomHugeTail(double obs_k, double n, td_real p_tdr, uint32_t complement, int32_t midp, uint32_t logp) {
  // obs_k and n are integers, n >= 2^52, 0 <= min(obs_k, n - obs_k) <= 2048
  // obs_k=n only possible when midp true
  // p or q could be extreme
  dd_real obs_k_ddr = ddr_maked(obs_k);
  dd_real p_ddr = ddr_make_td(p_tdr);
  dd_real q_ddr = ddr_negate(ddr_make_td(tdr_addd(p_tdr, -1.0)));
  if (complement) {
    obs_k_ddr = ddr_add2d(n, -obs_k);
    if (!midp) {
      // This may be lost to floating-point error.
      obs_k_ddr = ddr_addd(obs_k_ddr, -1);
    }
    swap_ddr(&p_ddr, &q_ddr);
  }
  double k = obs_k_ddr.x[0];
  if (ddr_is_zero(p_ddr)) {
    if ((k == 0) && midp) {
      return logp? -kLn2 : 0.5;
    }
    return logp? 0.0 : 1.0;
  }
  if (ddr_is_zero(q_ddr)) {
    if (k < n) {
      return logp? -INFINITY_D : 0.0;
    }
    assert(midp);
    return logp? -kLn2 : 0.5;
  }
  uint32_t calc_complement = 0;
  if (k > 2048) {
    // This updates obs_k_ddr/p_ddr/q_ddr.
    calc_complement = 1;
    if (!midp) {
      // Ensure that if complement and calc_complement are both true, the
      // previous -1 is lost to floating-point error iff this is.
      obs_k_ddr = ddr_addd(obs_k_ddr, 1);
    }
    k = ddr_negate(ddr_addd(obs_k_ddr, -n)).x[0];
    assert(k <= 2048);
    obs_k_ddr = ddr_maked(k);
    swap_ddr(&p_ddr, &q_ddr);
  }
  const double modal_k = floor(ddr_mul(p_ddr, ddr_add2d(n, 1)).x[0]);
  const double starting_k = MINV(modal_k, k);
  dd_real ln_prob_ddr = binom_ln_prob_loader(ddr_maked(starting_k), ddr_maked(n), p_ddr, q_ddr);
  dd_real lik_ddr = ddr_maked(1.0);
  dd_real sum_ddr;
  if (starting_k < k) {
    // We're actually to the right of the mode.  p must be tiny; need to be
    // careful about overflow/underflow.
    // Sum inwards from the mode.
    dd_real pdq_ddr = p_ddr;
    if (p_ddr.x[0] > 1.0 / (k2p800 * k2p100)) {
      // Don't want to worry about NaN here.
      pdq_ddr = ddr_accurate_div(p_ddr, q_ddr);
    }
    const double k_stop = k;
    k = starting_k;
    sum_ddr = ddr_maked(1);
    do {
      const dd_real nmk_ddr = ddr_add2d(n, -k);
      k += 1;
      lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_mul(pdq_ddr, nmk_ddr), k));
      if (k == k_stop) {
        if (midp) {
          lik_ddr = ddr_mul_pwr2(lik_ddr, 0.5);
        }
        sum_ddr = ddr_add(sum_ddr, lik_ddr);
        break;
      }
      sum_ddr = ddr_add(sum_ddr, lik_ddr);
    } while (lik_ddr.x[0] >= (k2m64 / 8));
    k = starting_k;
    lik_ddr = ddr_maked(1);
  } else {
    // Some obvious early-exit opportunities.
    if (((!logp) || calc_complement) && (ln_prob_ddr.x[0] < -1085 * kLn2)) {
      if (calc_complement) {
        return logp? 0.0 : 1.0;
      }
      return 0.0;
    }
    sum_ddr = ddr_maked(1.0 - 0.5 * midp);
  }
  // Now sum outward until next term is less than ~2^{-67} of starting term.

  // Next term is kq / ((n-k+1)p) times the current one.  Since k <= 2^11 and
  // n >= 2^52, q < 2^{-27} guarantees the first multiplier < 2^{-67}; worth
  // cheaply checking this up front so we don't have to worry about e.g.
  // qdp_ddr underflow.
  if ((k > 0) && (q_ddr.x[0] > 1.0 / (1 << 27))) {
    const dd_real qdp_ddr = ddr_accurate_div(q_ddr, p_ddr);
    do {
      lik_ddr = ddr_mul(lik_ddr, ddr_accurate_div(ddr_muld(qdp_ddr, k), ddr_add2d(n, -k)));
      k -= 1;
      sum_ddr = ddr_add(sum_ddr, lik_ddr);
    } while (lik_ddr.x[0] >= (k2m64 / 8));
  }
  if (!calc_complement) {
    ln_prob_ddr = ddr_add(ln_prob_ddr, ddr_log(sum_ddr));
    if (ln_prob_ddr.x[0] > 0) {
      ln_prob_ddr = ddr_maked(0);
    }
    if (logp) {
      return ln_prob_ddr.x[0];
    }
    return ddr_exp(ln_prob_ddr).x[0];
  }
  dd_real prob_ddr = ddr_negate(ddr_addd(ddr_mul(ddr_exp(ln_prob_ddr),
                                                 sum_ddr),
                                         -1));
  if (ddr_gtd(prob_ddr, 1.0)) {
    prob_ddr = ddr_maked(1);
  }
  if (!logp) {
    return prob_ddr.x[0];
  }
  return ddr_log(prob_ddr).x[0];
}

double PbinomExtremeSuccP(double obs_k, double n, td_real p_tdr, uint32_t complement, int32_t midp, uint32_t logp) {
  // min(p, 1-p) < 2^{-924}, n <= 2^900 (so n * min(p,1-p) < 2^{-24})
  // Need to be careful about underflow, but this case is otherwise
  // straightforward since the log-likelihood is either strictly and rapidly
  // decreasing or strictly and rapidly increasing.
  // If n >= 2^52, min(obs_k+1, n-obs_k) >= 40; this lets us simplify the
  // p<0.5 branch.
  dd_real p_ddr = ddr_make_td(p_tdr);
  dd_real q_ddr = ddr_make_td(tdr_negate(tdr_addd(p_tdr, -1.0)));
  dd_real k_ddr = ddr_maked(obs_k);
  dd_real nmk_ddr = ddr_add2d(n, -obs_k);
  if (complement) {
    k_ddr = nmk_ddr;
    nmk_ddr = ddr_maked(obs_k);
    if (!midp) {
      if (ddr_is_zero(k_ddr)) {
        return logp? -INFINITY_D : 0.0;
      }
      k_ddr = ddr_addd(k_ddr, -1);
      nmk_ddr = ddr_addd(nmk_ddr, 1);
    }
    swap_ddr(&p_ddr, &q_ddr);
  }
  if (p_ddr.x[0] < 0.5) {
    const double k = k_ddr.x[0];
    // pmf(0) = q^n
    // pmf(1) = p * q^{n-1} * n ~= np when q>(1 - 2^{-924}), n<2^52
    // pmf(2 or greater) = underflow unless n huge
    if (midp && (k == 0)) {
      // Since we can only get here when n < 2^52, difference from 0.5 is too
      // small to matter.
      return logp? -kLn2 : 0.5;
    }
    if (!logp) {
      // Difference from 1 only representable if we're in log-space.
      return 1.0;
    }
    // log(1-x) = -x - x^2/2 - ... ~= -x for small x.  If x underflows, final
    // result is indistinguishable from 1 even in log-space.
    if (k + (!midp) == 1) {
      // n < 2^52, p < 2^{-924}
      // log(cdf(0)) ~= -ccdf(0) ~= -pmf(1)
      // log(cdf(1) - 0.5 * pmf(1)) ~= -0.5 * pmf(1)
      return (0.5 * midp - 1) * ddr_muld(p_ddr, n).x[0];
    }
    // Either np < 2^{-24} and k >= 39, or np < 2^{-872} and k+(!midp) >= 2.
    // In both cases, we underflow.
    return 0;
  }

  if (ddr_is_zero(nmk_ddr)) {
    if (!midp) {
      return logp? 0 : 1;
    }
    return logp? -kLn2 : 0.5;
  }
  if (ddr_is_zero(q_ddr)) {
    return logp? -INFINITY_D : 0.0;
  }
  if (nmk_ddr.x[0] == 1) {
    const double pval = (1 - 0.5 * midp) * ddr_muld(q_ddr, n).x[0];
    return logp? log(pval) : pval;
  }
  // cdf(obs_k) guaranteed to underflow...
  if (!logp) {
    return 0;
  }
  // but log(cdf(obs_k)) does not, instead it's a large-magnitude negative
  // number.
  // log(cdf(obs_k)) = log(pmf(obs_k) + pmf(obs_k - 1) + ...)
  //                 = log(pmf(obs_k)) + log(1 + pmf(obs_k - 1)/pmf(obs_k) + ...)
  // pmf(obs_k) = p^{obs_k} q^{nmk} (n choose nmk)
  // pmf(obs_k - 1)/pmf(obs_k) = (q/p)(obs_k / (n - obs_k + 1))
  //                          ~= q(obs_k / (n - obs_k + 1))
  //                          <= 2^{-24} / 40, so this term matters in edge
  //                             case but subsequent terms too small
  dd_real retval_ddr = binom_ln_prob_loader(k_ddr, nmk_ddr, p_ddr, q_ddr);
  if (midp) {
    retval_ddr = ddr_sub(retval_ddr, _ddr_log2);
  }
  if (retval_ddr.x[0] > -(1LL << 38)) {
    // Pbinom() currently targets epsilon=2^{-67}.  If retval < -(2^38), this
    // term always has relative contribution smaller than that.  If not, it's
    // safe to compute this term with float64 precision.
    retval_ddr = ddr_addd(retval_ddr, q_ddr.x[0] * k_ddr.x[0] / (nmk_ddr.x[0] + 1));
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
    return binom_ln_prob_loader(ddr_maked(k), ddr_maked(n), p_ddr, q_ddr);
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

dd_real binom_ltail_lik_bfrac_ddr(double obs_k, double n, dd_real p_ddr, dd_real q_ddr) {
  double aa = obs_k + 1;
  double bb = n - obs_k;
  const dd_real orig_p_nmk_ddr = ddr_muld(p_ddr, bb);
  dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), orig_p_nmk_ddr);
  if (ay_minus_bx_ddr.x[0] < 0.0) {
    swap_f64(&aa, &bb);
    swap_ddr(&p_ddr, &q_ddr);
    ay_minus_bx_ddr = ddr_negate(ay_minus_bx_ddr);
  }
  return ddr_accurate_div(orig_p_nmk_ddr, ibeta_continued_fraction_ddr(ddr_maked(aa), ddr_maked(bb), n, p_ddr, q_ddr, ay_minus_bx_ddr));
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
