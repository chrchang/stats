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

#include "binom.h"

#include <assert.h>
#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Assumes 0 <= n < 2^52, k in [0, n].
double LnBinomCoeff(int64_t n, int64_t k) {
  if ((k == 0) || (k == n)) {
    return 0;
  }
  if (n < (1LL << 37)) {
    return ddr_sub(ddr_lfact(n),
                   ddr_add_lfacts(k, n-k)).x[0];
  }
  return tdr_sub(tdr_lfact(n),
                 tdr_add(tdr_lfact(k), tdr_lfact(n-k))).x[0];
}

// Ok to draw this line anywhere <= 2^39 (see error analysis below).
// I've set this to 2^36 since that's roughly where I could no longer easily
// find 1 ULP deviations from MPFR.
static inline uint32_t use_tdr_for_binom_lnprob(int64_t obs_tot) {
  return (obs_tot >= (1LL << 36));
}

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
//
// This implementation is a bit slow, but it's relatively simple and reliable.
// (ibeta_power_terms_d_ln() trades off a tiny bit of accuracy for a
// significant speed improvement.)
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

// Assumes 0 <= k <= n < 2^52, 0 < p < 1.
double BinomMass(int64_t k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t logp) {
  if (k == n) {
    k = 0;
    swap_ddr(&p_ddr, &q_ddr);
  }
  const dd_real ln_prob_ddr = binom_ln_prob_internal(k, n, p_ddr, q_ddr);
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  // Note that if we want to deliver j-bit precision for 0 < x < 2^{-512},
  // log(x) needs to be represented to ~(j+10)-bit precision.  (This
  // requirement has discouraged preexisting implementations from supporting
  // the logp and !logp cases with a shared core.)
  return ddr_exp(ln_prob_ddr).x[0];
}

// - succ_odds_ratio_tdr must be p/(1-p), where p is the expected success rate.
//
// - starting_lnprobv_ddr is expected to either be initialized to
//     log(succ_odds_ratio^obs_succ / (obs_succ! (obs_tot - obs_succ)!)),
//   or have x[0] initialized to DBL_MAX to indicate that that calculation
//   hasn't happened.  In the latter case, it may be set to the former value if
//   that is needed in the calculation.
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

// ibeta_fraction2_...() adapted from Boost 1.91.0.  This derived code is
// subject to the following license:
//
// *****
// Boost Software License - Version 1.0 - August 17th, 2003
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// *****

// static const double kLentzFpmin = DBL_MIN * 16;

static const double kLanczosDoubleSumDenom[13] = {0, 39916800, 120543840, 150917976, 105258076, 45995730, 13339535, 2637558, 357423, 32670, 1925, 66, 1};
static const double kLanczosDoubleSumExpgNumer[13] = {
  56906521.91347156388090791033559122686859,
  103794043.1163445451906271053616070238554,
  86363131.28813859145546927288977868422342,
  43338889.32467613834773723740590533316085,
  14605578.08768506808414169982791359218571,
  3481712.15498064590882071018964774556468,
  601859.6171681098786670226533699352302507,
  75999.29304014542649875303443598909137092,
  6955.999602515376140356310115515198987526,
  449.9445569063168119446858607650988409623,
  19.51992788247617482847860966235652136208,
  0.5098416655656676188125178644804694509993,
  0.006061842346248906525783753964555936883222
};

// this depends on the polynomial coefficients above
// exactly 808618867 * 2^{-27}, don't need to represent this as dd_real
static const double kLanczosDoubleG = 6.024680040776729583740234375;

double lanczos_sum_d_expg_scaled_imp(double zz, double* s2_ptr) {
  // zz currently guaranteed to be >1.
  zz = 1 / zz;
  double s1 = kLanczosDoubleSumExpgNumer[0];
  double s2 = kLanczosDoubleSumDenom[0];
  for (uint32_t uii = 1; uii != 13; ++uii) {
    s1 = prefer_fma(s1, zz, kLanczosDoubleSumExpgNumer[uii]);
    s2 = prefer_fma(s2, zz, kLanczosDoubleSumDenom[uii]);
  }
  *s2_ptr = s2;
  return s1;
}

dd_real ibeta_power_terms_d_ln(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real ay_minus_bx_ddr) {
  // Returns log((x^a)(y^b) / Beta(a,b))
  //       = log((x^a)(y^b)(a+b-1)! / ((a-1)!(b-1)!)).
  //
  // normalized always true
  // prefix always 1
  // aa and bb always large
  //
  // This actually holds up better than the obvious ddr-based implementation
  // when aa+bb approaches 2^52.
  double cc = aa + bb;
  const double gh = kLanczosDoubleG - 0.5;
  const dd_real agh_ddr = ddr_add2d(gh, aa);
  const dd_real bgh_ddr = ddr_add2d(gh, bb);
  const dd_real cgh_ddr = ddr_add2d(gh, cc);

  double numer_a;
  const double denom_a = lanczos_sum_d_expg_scaled_imp(aa, &numer_a);
  double numer_b;
  const double denom_b = lanczos_sum_d_expg_scaled_imp(bb, &numer_b);
  double denom_c;
  const double numer_c = lanczos_sum_d_expg_scaled_imp(cc, &denom_c);
  // Tried performing some of the following computations with float64s instead
  // of dd_reals, but that resulted in noticeably higher error than scipy.
  const dd_real term1_ddr = ddr_sqr(ddr_accurate_div(ddr_muld(ddr_mul2d(numer_a, numer_b), numer_c), ddr_muld(ddr_mul2d(denom_a, denom_b), denom_c)));
  const dd_real term2_ddr = ddr_accurate_div(ddr_mul(agh_ddr, bgh_ddr), ddr_mul(cgh_ddr, _ddr_e));
  dd_real result_ddr = ddr_mul_pwr2(ddr_log(ddr_mul(term1_ddr, term2_ddr)), 0.5);

  // Calculate l1 and l2 with extra precision, since magnitude can greatly
  // exceed that of ln(nonlog).
  // This removes the need for special cases.
  const dd_real l1_ddr = ddr_accurate_div(ddr_negate(ddr_add(ay_minus_bx_ddr, ddr_muld(q_ddr, gh))), agh_ddr);
  const dd_real l2_ddr = ddr_accurate_div(ddr_sub(ay_minus_bx_ddr, ddr_muld(p_ddr, gh)), bgh_ddr);
  return ddr_sort_and_add3(result_ddr, ddr_muld(ddr_log1p(l1_ddr), aa), ddr_muld(ddr_log1p(l2_ddr), bb));
}

double ibeta_continued_fraction_recip_d(double aa, double bb, double xx, double yy, dd_real ay_minus_bx_ddr, uint32_t inv, uint32_t midp_complement) {
  // see Boost continued_fraction_b()
  const double ay_minus_bx_plus1 = ay_minus_bx_ddr.x[0] + 1.0;
  double cc = (aa / (aa + 1.0)) * ay_minus_bx_plus1;
  const double two_minus_x = 2 - xx;
  // This provides a noticeable accuracy boost at reasonable computational
  // cost.
  dd_real ff_ddr = ddr_maked(cc);
  double dd = 0.0;
  double mm = 1.0;
  while (1) {
    const double denom = aa + 2 * mm - 1;
    // if xx is very small, precomputed xx * xx may underflow when actual
    // product here does not
    // (also possible for actual product to underflow)
    const double shared_middle_term = mm * (bb - mm) * xx / denom;
    const double cur_a = ((aa + mm - 1) / denom) * (aa + bb + mm - 1) * shared_middle_term * xx;
    double cur_b = prefer_fma((aa + mm) / (aa + 2 * mm + 1), prefer_fma(mm, two_minus_x, ay_minus_bx_plus1), mm + shared_middle_term);
    // cur_b += ((aa + mm) / (aa + 2 * mm + 1)) * prefer_fma(mm, two_minus_x, ay_minus_bx_plus1);
    mm += 1.0;
    dd = prefer_fma(cur_a, dd, cur_b);
    // Algorithm should terminate when cur_a decreases to 0 due to bb == mm or
    // underflow.  At and before that point, cur_b is always positive.
    /*
    if (dd == 0.0) {
      dd = kLentzFpmin;
    }
    */
    cc = cur_b + cur_a / cc;
    /*
    if (cc == 0.0) {
      cc = kLentzFpmin;
    }
    */
    dd = 1.0 / dd;
    const dd_real delta_ddr = ddr_mul2d(cc, dd);
    if (fabs(delta_ddr.x[0] - 1) <= k2m52) {
      double inv_ff = 1.0 / ff_ddr.x[0];
      if (midp_complement) {
        // If complement=0, inv=1 (a<->b not flipped):
        //   result_ln
        // = a log x + b log y + log((a+b-1)!) - log((a-1)!) - log((b-1)!)
        // = (k+1) log p + (n-k) log q + log(n!) - log(k!) - log((n-k-1)!)
        //
        //   log(0.5 * pmf(k))
        // = log 0.5 + k log p + (n-k) log q + log(n!) - log(k!) - log((n-k)!)
        // = log 0.5 - log p - log(n-k) + result_ln
        // = -log (2 * p * (n-k))) + result_ln
        //
        //   log(exp(result_ln - log(ff)) - 0.5 * pmf(k))
        // = log(exp(result_ln)/ff - exp(result_ln - log(2p(n-k))))
        // = log(result/ff - result/(2p(n-k)))
        // = log(result * (1/ff - 0.5/(p(n-k))))
        // = result_ln + log(1/ff - 0.5/(p(n-k)))
        double signed_p_nmk;
        if (midp_complement == inv + 1) {
          // aa<->bb, xx<->yy flip was performed earlier.
          signed_p_nmk = -aa * yy;
        } else {
          signed_p_nmk = bb * xx;
        }
        inv_ff += 0.5 / signed_p_nmk;
      }
      return inv_ff;
    }
    ff_ddr = ddr_mul(ff_ddr, delta_ddr);
  }
}

// Adaptations of DiDonato and Morris's BFRAC, which is in turn based on a
// continued fraction introduced in
//   Aroian LA (1941) Continued fractions for the incomplete beta function.
//   Annals of Mathematical Statistics, 12.
// For most larger cases, this continued fraction converges more quickly than
// binomial partial sums.  The _ln_ddr1 function makes limited use of dd_real
// precision to address the worst precision bottlenecks, while the _ddr2
// function trades off speed for provably great precision.
//
// (I still have work to do in understanding the derivation and properties of
// this continued fraction well enough to take a real shot at improving e.g.
// the rather similar hypergeometric cdf calculation.)
dd_real ibeta_fraction2_ln_ddr1(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real ay_minus_bx_ddr, uint32_t inv, uint32_t midp_complement) {
  // normalized always true, min(aa, bb) >= 40, max much larger
  // (this should still yield correct results for smaller min(aa, bb), but it
  // looks relatively inefficient in that case.  todo: benchmark.)
  // caller responsible for guaranteeing ay - bx >= 0
  dd_real result_ln_ddr = ibeta_power_terms_d_ln(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr);
  const double result_incr = log(ibeta_continued_fraction_recip_d(aa, bb, p_ddr.x[0], q_ddr.x[0], ay_minus_bx_ddr, inv, midp_complement));
  result_ln_ddr = ddr_addd(result_ln_ddr, result_incr);
  if (!inv) {
    return result_ln_ddr;
  }
  return ddr_log1p(ddr_negate(ddr_exp(result_ln_ddr)));
}

dd_real ibeta_continued_fraction_ddr(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real ay_minus_bx_ddr) {
  const double cf_eps = k2m64 * (1.0 / (1 << 3));
  const dd_real ay_minus_bx_plus1_ddr = ddr_addd(ay_minus_bx_ddr, 1);
  dd_real cc_ddr = ddr_divd(ddr_muld(ay_minus_bx_plus1_ddr, aa), aa + 1);
  const dd_real two_minus_x_ddr = ddr_addd(q_ddr, 1);
  dd_real ff_ddr = cc_ddr;
  dd_real dd_ddr = ddr_maked(0);
  double mm = 1.0;
  while (1) {
    const double denom = aa + 2 * mm - 1;
    const dd_real shared_middle_term_ddr = ddr_divd(ddr_mul(ddr_mul2d(mm, bb - mm), p_ddr), denom);
    // if p_ddr is very small, precomputed p_ddr * p_ddr may underflow when
    // actual product here does not
    // (also possible for actual product to underflow)
    // (we probably just want to handle tiny p in its own function?)
    const dd_real cur_a_ddr = ddr_divd(ddr_mul(ddr_mul(ddr_mul2d(aa + mm - 1, aa + bb + mm - 1), shared_middle_term_ddr), p_ddr), denom);
    dd_real cur_b_ddr = ddr_addd(shared_middle_term_ddr, mm);
    cur_b_ddr = ddr_add(cur_b_ddr, ddr_divd(ddr_muld(ddr_add(ddr_muld(two_minus_x_ddr, mm), ay_minus_bx_plus1_ddr), aa + mm), aa + 2 * mm + 1));
    mm += 1.0;
    dd_ddr = ddr_add(ddr_mul(cur_a_ddr, dd_ddr), cur_b_ddr);
    // Algorithm should terminate when cur_a decreases to 0 due to bb == mm or
    // underflow.  At and before that point, cur_b is always positive.
    /*
    if (dd == 0.0) {
      dd = kLentzFpmin;
    }
    */
    cc_ddr = ddr_add(cur_b_ddr, ddr_accurate_div(cur_a_ddr, cc_ddr));
    /*
    if (cc == 0.0) {
      cc = kLentzFpmin;
    }
    */
    dd_ddr = ddr_accurate_div(ddr_maked(1.0), dd_ddr);
    const dd_real delta_ddr = ddr_mul(cc_ddr, dd_ddr);
    // If I correctly understand what's going on here, (delta - 1) has
    // alternating sign and decreasing magnitude, so this should ensure less
    // than cf_eps relative error is coming from incomplete evaluation of the
    // continued fraction.  (Recall that we actually need to limit the
    // logarithm's absolute error to < ~2^{-63} to achieve the desired level of
    // accuracy when DBL_MIN < p < e^{-512} and logp=False.)
    //
    // Worst case, this takes around a million iterations.
    //
    // Update (3 Jun 2026): cf_eps tightened to 2^{-67} after comparing against
    // MPFR.  We want 1 ULP errors to require a non-negligible amount of work
    // to find.  The resulting ~5-25% speed hit is an acceptable price to pay
    // for a still-much-faster-than-MPFR function that can be treated as
    // baseline truth for most other purposes.
    //
    // The use_tdr_for_binom_lnprob() threshold should also be lowered if we
    // want to push cf_eps below 2^{-67}.
    if ((delta_ddr.x[0] == 1.0) && (fabs(delta_ddr.x[1]) <= cf_eps)) {
      return ff_ddr;
    }
    ff_ddr = ddr_mul(ff_ddr, delta_ddr);
  }
}

double ibeta_fraction2_ddr2(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real ay_minus_bx_ddr, uint32_t inv, uint32_t logp) {
  // (a is a synonym for k+1, b is a synonym for n-k, x is a synonym for p, y
  // is a synonym for q)
  //   log((x^a)(y^b) / Beta(a,b))
  // = a log x + b log y + log((a+b-1)!) - log((a-1)!) - log((b-1)!)
  const double a_plus_b = aa + bb;
  const uint32_t p_is_half = ddr_is(p_ddr, 0.5);
  dd_real result_ln_ddr;
  if (!use_tdr_for_binom_lnprob(S_CAST(int64_t, a_plus_b - 1))) {
    dd_real ddrs[5];
    ddrs[0] = ddr_lfact(a_plus_b - 1);
    ddrs[1] = ddr_negate(ddr_lfact(aa - 1));
    ddrs[2] = ddr_negate(ddr_lfact(bb - 1));
    if (p_is_half) {
      ddrs[3] = ddr_muld(_ddr_log05, a_plus_b);
    } else {
      ddrs[3] = ddr_muld(ddr_log_2arg(p_ddr, q_ddr), aa);
      ddrs[4] = ddr_muld(ddr_log_2arg(q_ddr, p_ddr), bb);
    }
    result_ln_ddr = ddr_sort_and_add(5 - p_is_half, ddrs);
  } else {
    td_real tdrs[5];
    tdrs[0] = tdr_lfact(a_plus_b - 1);
    tdrs[1] = tdr_negate(tdr_lfact(aa - 1));
    tdrs[2] = tdr_negate(tdr_lfact(bb - 1));
    if (p_is_half) {
      tdrs[3] = tdr_muld(_tdr_log05, a_plus_b);
    } else {
      if (ddr_ltd(p_ddr, 0.5)) {
        tdrs[3] = tdr_muld(tdr_log(tdr_make_dd(p_ddr)), aa);
        tdrs[4] = tdr_muld(tdr_log1p(tdr_make_dd(ddr_negate(p_ddr))), bb);
      } else {
        tdrs[3] = tdr_muld(tdr_log1p(tdr_make_dd(ddr_negate(q_ddr))), aa);
        tdrs[4] = tdr_muld(tdr_log(tdr_make_dd(q_ddr)), bb);
      }
    }
    result_ln_ddr = ddr_make_td(tdr_sort_and_add(5 - p_is_half, tdrs));
  }
  // Could tighten this bound.
  if ((result_ln_ddr.x[0] < -1418.0) && (inv || (!logp))) {
    return (logp || (!inv))? 0.0 : 1.0;
  }
  const dd_real ff_ddr = ibeta_continued_fraction_ddr(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr);
  result_ln_ddr = ddr_sub(result_ln_ddr, ddr_log(ff_ddr));
  if (!inv) {
    return logp? result_ln_ddr.x[0] : ddr_exp(result_ln_ddr).x[0];
  }
  const dd_real neg_result_ddr = ddr_negate(ddr_exp(result_ln_ddr));
  return logp? ddr_log1p(neg_result_ddr).x[0] : ddr_addd(neg_result_ddr, 1.0).x[0];
}

// Requires 0 <= obs_succ <= obs_tot < 2^52 and 2^{-960} < p < 1.
// (Sometimes works for 0 < p <= 2^{-960}, but let's leave that out of the
// function contract until the holes in that region are plugged in.)
//
// Benchmark results revealed that Boost 1.91 ibetac(k+1, n-k, p) (which is
// called by scipy.stats.binom.logcdf()) became faster than this function's
// initial implementation once obs_tot was ~1000, and its results were
// acceptably accurate.  We now use its main algorithm when obs_tot > 512 and
// min(k+1, n-k) >= 40.
//
// Interestingly, the scipy implementation has much higher overhead, even after
// initialization, despite relying on very similar C++ code.  E.g.
//
//   >>> import exact_tests, scipy, timeit
//   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
//   0.01796633400954306
//   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
//   0.0219896660419181
//   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
//   0.017003458982799202
//   >>> timeit.timeit(lambda: scipy.stats.binom.logcdf(157000000, 419430500, 0.375), number=10000)
//   1.005605333019048
//   >>> timeit.timeit(lambda: scipy.stats.binom.logcdf(157000000, 419430500, 0.375), number=10000)
//   0.22769695916213095
//   >>> timeit.timeit(lambda: scipy.stats.binom.logcdf(157000000, 419430500, 0.375), number=10000)
//   0.24197920807637274
//   >>> timeit.timeit(lambda: scipy.stats.binom.cdf(157000000, 419430500, 0.375), number=10000)
//   0.27701791608706117
//   >>> timeit.timeit(lambda: scipy.stats.binom.cdf(157000000, 419430500, 0.375), number=10000)
//   0.27656474988907576
//   >>> timeit.timeit(lambda: scipy.stats.binom.cdf(157000000, 419430500, 0.375), number=10000)
//   0.27791083394549787

// See Pbinom() below for a higher-accuracy variant of this function; this one
// is limited by the float64 precision of the Lanczos and
// continued-fraction-coefficient calculations, and everything else here is
// tuned to that level of relative error.
double PbinomApprox(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t complement, int32_t midp, uint32_t logp) {
  if ((obs_k < 0) || (obs_k > n)) {
    if ((obs_k < 0) == complement) {
      return logp? 0.0 : 1.0;
    }
    return logp? (0.0 / 0.0) : 0.0;
  }
  if ((n > 512) && (MINV(obs_k + 1, n - obs_k) >= 40)) {
    double aa = obs_k + 1;
    double bb = n - obs_k;
    dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), ddr_muld(p_ddr, bb));
    uint32_t inv = !complement;
    if (ay_minus_bx_ddr.x[0] < 0.0) {
      swap_f64(&aa, &bb);
      swap_ddr(&p_ddr, &q_ddr);
      ay_minus_bx_ddr = ddr_negate(ay_minus_bx_ddr);
      inv = !inv;
    }
    // (took a brief look at Boost's use_asym branch, don't see a reasonable
    // way to use it without sacrificing accuracy, and current code seems fast
    // enough.)
    dd_real result_ln_ddr = ibeta_fraction2_ln_ddr1(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr, inv, midp * (1 + complement));
    if (logp) {
      return result_ln_ddr.x[0];
    }
    return ddr_exp(result_ln_ddr).x[0];
  }
  if (complement) {
    obs_k = n - obs_k - (!midp);
    if (obs_k < 0) {
      return logp? (0.0 / 0.0) : 0.0;
    }
    swap_ddr(&p_ddr, &q_ddr);
  }
  const double pdq = ddr_accurate_div(p_ddr, q_ddr).x[0];
  double k = obs_k;
  double nmk = n - obs_k;
  if (k > nmk * pdq) {
    // We're at or to the right of the mode.
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know the p-value > 1 - logp? DBL_MIN : 2^{-54}
    // (at which point we just return log(1) or 1; in the logp case, don't want
    // to risk imposing a surprising denormal-handling performance penalty for
    // no good reason), or remaining left likelihoods are smaller than the
    // precision limit.
    const double first_right_mult = pdq * nmk / (k + 1);
    // r + r^2 + ... = r / (1-r)
    const double right_upper_bound = 0.5 * midp + first_right_mult / (1 - first_right_mult);
    if (right_upper_bound == 0.0) {
      // p-value is exactly 1 when nmk==0 and midp is false
      return logp? 0 : 1;
    }

    // Scale our starting likelihood so that we overflow to INFINITY when we'd
    // want to early-exit and return log(1) or 1; this saves us a comparison in
    // the loop.
    const double start_lik = (DBL_MAX * (logp? DBL_MIN : k2m54)) / right_upper_bound;
    double lik = start_lik;
    double left_sum = start_lik;
    while (1) {
      nmk += 1;
      lik *= k / (pdq * nmk);
      k -= 1;
      const double preadd = left_sum;
      left_sum += lik;
      if (left_sum == preadd) {
        break;
      }
    }
    if (left_sum == INFINITY_D) {
      return logp? 0 : 1;
    }
    // Now compute the right-sum to the precision limit.
    double right_sum = first_right_mult * start_lik;
    k = obs_k + 1;
    nmk = n - obs_k - 1;
    lik = right_sum;
    while (1) {
      k += 1;
      lik *= pdq * nmk / k;
      nmk -= 1;
      const double preadd = right_sum;
      right_sum += lik;
      if (right_sum == preadd) {
        break;
      }
    }
    // For one-sided test, slightly more convenient to exclude midp term from
    // left_sum and right_sum since it just cancels out in denom
    const double midp_numer = -0.5 * midp * start_lik;
    const double denom = right_sum + left_sum;
    if (!logp) {
      return (left_sum + midp_numer) / denom;
    }
    return log1p((midp_numer - right_sum) / denom);
  }
  // We're to the left of the mode, and are responsible for tiny p-values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise, we evaluate the starting log-likelihood with dd_reals to work
  // around catastrophic cancellation, and then iterate leftward to the
  // precision limit.
  const double ln_first_inward_mult = log(pdq * nmk / (k + 1));
  double overflow_steps_lower_bound = 1LL << 52;
  // log(DBL_MAX / 2^52) = 673.739...
  if (ln_first_inward_mult > 673.739 * k2m52) {
    overflow_steps_lower_bound = 673.739 / ln_first_inward_mult;
  }
  const double nd = n;
  const double modal_k = nd * pdq / (1 + pdq);
  if (k + overflow_steps_lower_bound > modal_k) {
    double lik = 1;
    double right_sum = 0;
    while (1) {
      k += 1;
      lik *= pdq * nmk / k;
      nmk -= 1;
      const double preadd = right_sum;
      right_sum += lik;
      if (right_sum == preadd) {
        break;
      }
    }
    k = obs_k;
    nmk = n - obs_k;
    lik = 1;
    double left_sum = 1;
    while (1) {
      nmk += 1;
      lik *= k / (pdq * nmk);
      k -= 1;
      const double preadd = left_sum;
      left_sum += lik;
      if (left_sum == preadd) {
        break;
      }
    }
    const double pval = (left_sum - 0.5 * midp) / (left_sum + right_sum);
    return logp? log(pval) : pval;
  }
  const dd_real pdq_ddr = ddr_maked(pdq);
  const dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(ddr_log(pdq_ddr), k),
            ddr_add_lfacts(k, nmk));
  dd_real ln_nmk_ddr = _ddr_log05;
  if (pdq != 1.0) {
    // log(1 / (1 + pdq)) = -log(1 + pdq)
    ln_nmk_ddr = ddr_negate(ddr_log(ddr_addd(pdq_ddr, 1.0)));
  }
  const dd_real lnprobf_ddr =
    ddr_add(ddr_lfact(nd), ddr_muld(ln_nmk_ddr, nd));
  const dd_real starting_lnprob_ddr = ddr_add(lnprobf_ddr, starting_lnprobv_ddr);
  // left_sum is the sum of < 2^52 terms, each of which is <= 1, so if
  // starting_lnprob < DBL_MIN / 2^52, final return value should always be 0
  // when logp=false and we're flushing denormals to zero.  DBL_MIN is
  // 2^{-1022}.
  //
  // 2^{-1074} is the smallest positive denormal, and (1 + epsilon) * 2^{-1075}
  // is the smallest number that should be rounded up to it, so -1074 can be
  // replaced with -1127 if we want this function to return denormals.
  //
  // (Yes, a tighter bound could be established for left_sum if it matters.)
  if ((!logp) && (starting_lnprob_ddr.x[0] < -1074 * kLn2)) {
    return 0;
  }
  double lik = 1;
  double left_sum = 1 - 0.5 * midp;
  while (1) {
    nmk += 1;
    lik *= k / (pdq * nmk);
    k -= 1;
    const double preadd = left_sum;
    left_sum += lik;
    if (left_sum == preadd) {
      break;
    }
  }
  return join_log_and_nonlog(starting_lnprob_ddr, left_sum, logp);
}

// Assumes 0 <= n < 2^52, 2^{-960} <= p < 1.
// Should consistently achieve <1 ULP relative error; and practically always
// <0.6 ULP except when n is close to 2^52.
//
// See PbinomApprox() above for the faster variant of this function which
// doesn't try to get the last few bits right.
double Pbinom(int64_t obs_k, int64_t n, dd_real p_ddr, dd_real q_ddr, uint32_t complement, uint32_t logp) {
  if ((obs_k < 0) || (obs_k > n - S_CAST(int64_t, complement))) {
    if ((obs_k < 0) == complement) {
      return logp? 0.0 : 1.0;
    }
    return logp? (0.0 / 0.0) : 0.0;
  }
  // Benchmarked various values of both thresholds, this seems good on my Mac
  if ((n > 131072) && (MINV(obs_k, n - obs_k) >= 2048)) {
    double aa = obs_k + 1;
    double bb = n - obs_k;
    dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), ddr_muld(p_ddr, bb));
    uint32_t inv = !complement;
    if (ay_minus_bx_ddr.x[0] < 0.0) {
      swap_f64(&aa, &bb);
      swap_ddr(&p_ddr, &q_ddr);
      ay_minus_bx_ddr = ddr_negate(ay_minus_bx_ddr);
      inv = !inv;
    }
    return ibeta_fraction2_ddr2(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr, inv, logp);
  }
  if (complement) {
    obs_k = n - obs_k - 1;
    swap_ddr(&p_ddr, &q_ddr);
  }
  const dd_real pdq_ddr = ddr_accurate_div(p_ddr, q_ddr);
  const dd_real qdp_ddr = ddr_accurate_div(q_ddr, p_ddr);
  double k = obs_k;
  double nmk = n - obs_k;
  if (k > nmk * pdq_ddr.x[0]) {
    // We're at or to the right of the mode.  Don't have to worry about cmf
    // values < DBL_MIN, so even for e.g. n > 2^51 our relative error is not
    // inflated by a need to work in log-space.

    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know we can safely return a cmf value of 1, or
    // remaining left likelihoods are smaller than the precision limit.
    const dd_real first_right_mult_ddr = ddr_divd(ddr_muld(pdq_ddr, nmk), k+1);
    const double right_upper_bound = first_right_mult_ddr.x[0] / ddr_negate(ddr_subd(first_right_mult_ddr, 1.0)).x[0];
    // Function contract ensures right_upper_bound > 0.

    // Scale our starting likelihood so that we overflow to INFINITY or nan
    // when we'd want to early-exit and return 1; this saves us a comparison in
    // the loop.
    const double start_lik = (DBL_MAX * (logp? DBL_MIN : k2m54)) / right_upper_bound;
    // We want to compute left_sum_ddr to at most 2^{-67} relative error (see
    // cf_eps discussion in ibeta_fraction_ddr2()), but we also want to drop
    // down from this slow dd_real-based loop to the much faster float64-based
    // loop as soon as we can prove that won't make us miss the accuracy
    // target.
    //
    // The lastp values we add to left_tail_sum have error bounded above by 2.5
    // ULPs, 4.5 ULPs, 6.5 ULPs, etc.  Thus, left_tail_sum should have error
    // bounded above by 3 ULPs, then 8, 15, 24, 35, ...; and k^2 is a loose
    // upper bound on the final error.
    //
    // (todo: actually, we can bound with 2.01, 3.51, 5.01, 6.51, 8.01, ..., so
    // we should be able to throw a ~3/4 multiplier in front?)
    //
    // 2^{-67} corresponds to at least 2^{52-67} ULPs.
    const double min_incr_left = (1.0 / (1 << 15)) / (k * k);
    dd_real lik_ddr = ddr_maked(start_lik);
    dd_real left_sum_ddr = lik_ddr;
    do {
      nmk += 1;
      lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
      k -= 1;
      left_sum_ddr = ddr_add(left_sum_ddr, lik_ddr);
      // This overflows to nan rather than INFINITY in my testing; I'll write
      // this to be agnostic to that detail.
    } while (lik_ddr.x[0] > left_sum_ddr.x[0] * min_incr_left);
    if (!(left_sum_ddr.x[0] < INFINITY_D)) {
      return logp? 0.0 : 1.0;
    }
    if (k > 0) {
      // Continue the calculation with ordinary precision.
      const double qdp = qdp_ddr.x[0];
      double lik = lik_ddr.x[0];
      double left_tail_sum = 0.0;
      while (1) {
        nmk += 1;
        lik *= qdp * k / nmk;
        k -= 1;
        const double preadd = left_tail_sum;
        left_tail_sum += lik;
        if (left_tail_sum == preadd) {
          break;
        }
      }
      left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
      if (!(left_sum_ddr.x[0] < INFINITY_D)) {
        return logp? 0.0 : 1.0;
      }
    }

    // Now compute the right-sum.
    dd_real right_sum_ddr = ddr_muld(first_right_mult_ddr, start_lik);
    k = obs_k + 1;
    nmk = n - obs_k - 1;
    lik_ddr = right_sum_ddr;
    if (nmk > 0) {
      const double min_incr_right = (1.0 / (1 << 15)) / (nmk * nmk);
      do {
        k += 1;
        lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(pdq_ddr, nmk), k));
        nmk -= 1;
        right_sum_ddr = ddr_add(right_sum_ddr, lik_ddr);
      } while (lik_ddr.x[0] > right_sum_ddr.x[0] * min_incr_right);
      if (nmk > 0) {
        // Continue the calculation with ordinary precision.
        const double pdq = pdq_ddr.x[0];
        double lik = lik_ddr.x[0];
        double right_tail_sum = 0.0;
        while (1) {
          k += 1;
          lik *= pdq * nmk / k;
          nmk -= 1;
          const double preadd = right_tail_sum;
          right_tail_sum += lik;
          if (right_tail_sum == preadd) {
            break;
          }
        }
        right_sum_ddr = ddr_addd(right_sum_ddr, right_tail_sum);
      }
    }
    const dd_real denom_ddr = ddr_add(left_sum_ddr, right_sum_ddr);
    if (!(denom_ddr.x[0] < INFINITY_D)) {
      return logp? 0.0 : 1.0;
    }
    const dd_real one_minus_prob_ddr = ddr_accurate_div(right_sum_ddr, denom_ddr);
    if (!logp) {
      return -ddr_subd(one_minus_prob_ddr, 1.0).x[0];
    }
    return ddr_log1p(ddr_negate(one_minus_prob_ddr)).x[0];
  }
  // We're at or to the left of the mode, and are responsible for tiny cmf
  // values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise, we use ddr_lfact() and friends to compute the starting
  // log-likelihood (this is the only step that could introduce >0.1 ULP error,
  // and then only for huge n), and accumulate the tail-sum from there.
  const double ln_first_inward_mult = log(pdq_ddr.x[0] * nmk / (k + 1));
  double overflow_steps_lower_bound = 1LL << 52;
  // log(DBL_MAX / 2^52) = 673.739...
  if (ln_first_inward_mult > 673.739 * k2m52) {
    overflow_steps_lower_bound = 673.739 / ln_first_inward_mult;
  }
  const double nd = n;
  const double modal_k = nd * p_ddr.x[0];
  if (k + overflow_steps_lower_bound > modal_k) {
    // (ok, this duplicated code belongs in its own function...)
    dd_real lik_ddr = ddr_maked(1.0);
    dd_real right_sum_ddr = ddr_maked(0.0);
    if (nmk > 0) {
      // Note that the k=28, n=56, p=0.5, logp=False case can ~randomly
      // mismatch MPFR by 1 ULP when this value is perturbed, because the true
      // cdf value is exactly on the boundary between two float64s.  (And on
      // such a mismatch, it's possible that it's the MPFR implementation
      // that isn't rounding toward even, we're lucky that it rounds correctly
      // in this case.)
      //
      // If others start treating this function as a source of ground truth, it
      // may be reasonable to use uint64 arithmetic to handle cases where n and
      // p are such that all cdf values are multiples of 2^{-64}.  That should
      // nail most exactly-on-the-boundary cases which naturally arise in
      // practice.
      const double min_incr_right = (1.0 / (1 << 15)) / (nmk * nmk);
      do {
        k += 1;
        lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(pdq_ddr, nmk), k));
        nmk -= 1;
        right_sum_ddr = ddr_add(right_sum_ddr, lik_ddr);
      } while (lik_ddr.x[0] > right_sum_ddr.x[0] * min_incr_right);
      if (nmk > 0) {
        const double pdq = pdq_ddr.x[0];
        double lik = lik_ddr.x[0];
        double right_tail_sum = 0.0;
        while (1) {
          k += 1;
          lik *= pdq * nmk / k;
          nmk -= 1;
          const double preadd = right_tail_sum;
          right_tail_sum += lik;
          if (right_tail_sum == preadd) {
            break;
          }
        }
        right_sum_ddr = ddr_addd(right_sum_ddr, right_tail_sum);
      }
    }
    k = obs_k;
    nmk = n - obs_k;
    lik_ddr = ddr_maked(1.0);
    dd_real left_sum_ddr = lik_ddr;
    if (k > 0) {
      const double min_incr_left = (1.0 / (1 << 15)) / (k * k);
      do {
        nmk += 1;
        lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
        k -= 1;
        left_sum_ddr = ddr_add(left_sum_ddr, lik_ddr);
      } while (lik_ddr.x[0] > left_sum_ddr.x[0] * min_incr_left);
      if (k > 0) {
        const double pdq = pdq_ddr.x[0];
        double lik = lik_ddr.x[0];
        double left_tail_sum = 0.0;
        while (1) {
          nmk += 1;
          lik *= k / (pdq * nmk);
          k -= 1;
          const double preadd = left_tail_sum;
          left_tail_sum += lik;
          if (left_tail_sum == preadd) {
            break;
          }
        }
        left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
      }
    }
    const dd_real prob_ddr = ddr_accurate_div(left_sum_ddr, ddr_add(left_sum_ddr, right_sum_ddr));
    if (!logp) {
      return prob_ddr.x[0];
    }
    return ddr_log(prob_ddr).x[0];
  }
  dd_real ln_prob_ddr = binom_ln_prob_internal(obs_k, n, p_ddr, q_ddr);
  if ((!logp) && (ln_prob_ddr.x[0] < -1074 * kLn2)) {
    return 0.0;
  }
  dd_real lik_ddr = ddr_maked(1.0);
  dd_real left_sum_ddr = lik_ddr;
  if (k > 0) {
    const double min_incr_left = (1.0 / (1 << 15)) / (k * k);
    do {
      nmk += 1;
      lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
      k -= 1;
      left_sum_ddr = ddr_add(left_sum_ddr, lik_ddr);
    } while (lik_ddr.x[0] > left_sum_ddr.x[0] * min_incr_left);
    if (k > 0) {
      // Continue the calculation with ordinary precision.
      const double qdp = qdp_ddr.x[0];
      double lik = lik_ddr.x[0];
      double left_tail_sum = 0.0;
      while (1) {
        nmk += 1;
        lik *= qdp * k / nmk;
        k -= 1;
        const double preadd = left_tail_sum;
        left_tail_sum += lik;
        if (left_tail_sum == preadd) {
          break;
        }
      }
      left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
    }
  }
  ln_prob_ddr = ddr_add(ln_prob_ddr, ddr_log(left_sum_ddr));
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  return ddr_exp(ln_prob_ddr).x[0];
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

dd_real binom_ltail_lik_simple_ddr(double k, double nmk, dd_real lik_ddr, dd_real qdp_ddr, double allowed_ulp_err) {
  dd_real tailsum_ddr = lik_ddr;
  if (k > 0) {
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
  }
  return tailsum_ddr;
}

// Returns smallest nonnegative k for which cdf(k) >= targetp.
// Assumes 0 <= n < 2^52, and should achieve targetp relative error < 2^{-54}.
// The goal is to support a higher-level qbinom() function where e.g.
//   qbinom(fractions.Fraction(1, 9), 2, fractions.Fraction(2, 3))
// can be trusted to be 0, and
//   qbinom(fractions.Fraction(2**53, (2**53 - 1) * 9), 2, fractions.Fraction(2, 3))
// can be trusted to be 1.
// The QbinomHalfUlp() entry point subtracts the natural epsilon value (0.5
// times the value of the least significant bit in targetp_or_lnp_ddr.x[0]) off
// of targetp_or_lnp_ddr, for the goal described above, before calling
// Qbinom().
//
// Probable todo: (except possibly on some very large cases) start with faster
// lower-accuracy interval-math calculation, and fall back on reliable
// high-accuracy calculation only when needed.
int64_t Qbinom(dd_real targetp_or_lnp_ddr, int64_t n, dd_real succp_ddr, uint32_t log_target) {
  if (ddr_is_zero(targetp_or_lnp_ddr) || ddr_is_zero(succp_ddr) || (n == 0)) {
    return 0;
  }
  if (ddr_is(succp_ddr, 1) || ddr_is(targetp_or_lnp_ddr, 1)) {
    return n;
  }
  // If targetp > 0.5, invert.
  const uint32_t inv = ((!log_target) && (targetp_or_lnp_ddr.x[0] > 0.5)) || (log_target && (targetp_or_lnp_ddr.x[0] > _ddr_log05.x[0]));
  dd_real failp_ddr = ddr_negate(ddr_subd(succp_ddr, 1.0));
  if (inv) {
    swap_ddr(&succp_ddr, &failp_ddr);
    if (!log_target) {
      targetp_or_lnp_ddr = ddr_negate(ddr_subd(targetp_or_lnp_ddr, 1.0));
    } else {
      targetp_or_lnp_ddr = ddr_negate(ddr_expm1(targetp_or_lnp_ddr));
      log_target = 0;
    }
  }
  // n=1 doesn't play well with the current initial-guess algorithm, and is
  // straightforward to handle directly.
  if (n == 1) {
    // cdf(0) = failp.
    const int64_t k = ddr_lt(failp_ddr, targetp_or_lnp_ddr);
    return inv? (1 - k) : k;
  }
  // We make an initial guess, use Newton's method to refine it (in the same
  // way as BinomTwoSidedP() when jumping from one tail to the other),
  // calculate tail probability to sufficient accuracy, and then start moving k
  // inward and updating tail probability until it crosses targetp.  This
  // avoids repetition of the tail-probability calculation.
  //
  // Initial guess is based on fitting a quadratic to three points of the
  // log-probability function near the mode.  Log-likelihood is
  //   log(p^k (1-p)^{n-k} n! / (k! (n-k)!))
  //   = k log(p) + (n-k) log(1-p) + log(n!) - log(k!) - log((n-k)!)
  //   = <constant> + k(log(p)-log(1-p)) - log(gamma(k+1)) - log(gamma(n+1-k))
  // First derivative w.r.t. k is
  //   log(p) - log(1-p) - digamma(k+1) + digamma(n+1-k)
  // Second derivative is
  //   -trigamma(k+1) - trigamma(n+1-k)
  //   ~= -(1/(k+1) + 1/(n+1-k))
  // which is slowly varying for large n and k.
  //
  // Probable todo: compare to R qbinom()'s use of Cornish-Fisher expansion.  I
  // expect the R approach to be better on this front.
  int64_t mode = S_CAST(int64_t, prev_float64((n + 1) * succp_ddr.x[0]));
  mode = mode + (mode == 0) - (mode == n);
  const uint32_t use_tdr = (n >= (1LL << 39));
  td_real n_lfact_tdr;
  td_real logp_tdr;
  td_real logq_tdr;
  if (!use_tdr) {
    n_lfact_tdr = tdr_make_dd(ddr_lfact(n));
    logp_tdr = tdr_make_dd(ddr_log(succp_ddr));
    logq_tdr = tdr_make_dd(ddr_log(failp_ddr));
  } else {
    n_lfact_tdr = tdr_lfact(n);
    logp_tdr = tdr_log(tdr_make_dd(succp_ddr));
    logq_tdr = tdr_log(tdr_make_dd(failp_ddr));
  }
  const dd_real pdq_ddr = ddr_accurate_div(succp_ddr, failp_ddr);
  const dd_real qdp_ddr = ddr_accurate_div(failp_ddr, succp_ddr);
  // This is usually overkill, but if n is huge we need to compute these to
  // high precision to make a decent initial guess.
  const dd_real mode_lnprob_ddr = ddr_sub(ddr_sort_and_add3(ddr_muld(ddr_make_td(logp_tdr), mode), ddr_muld(ddr_make_td(logq_tdr), n - mode), ddr_make_td(n_lfact_tdr)),
                                          ddr_add_lfacts(mode, n - mode));
  const dd_real modem1_lnprob_incr_ddr = ddr_log(ddr_divd(ddr_muld(qdp_ddr, mode), n + 1 - mode));
  const dd_real modep1_lnprob_incr_ddr = ddr_log(ddr_divd(ddr_muld(pdq_ddr, n - mode), mode + 1));

  // For better numerical stability, express this quadratic in terms of
  //   x := k - mode
  // instead of x := k.
  const double x2_coeff = 0.5 * ddr_add(modem1_lnprob_incr_ddr, modep1_lnprob_incr_ddr).x[0];
  const double x1_coeff = 0.5 * ddr_sub(modep1_lnprob_incr_ddr, modem1_lnprob_incr_ddr).x[0];
  const double x0_coeff = mode_lnprob_ddr.x[0];

  // 1. Identify k<mode such that
  //      pmf(k) * (k+1) < targetp
  //    If no such k>0 exists, initialize k=0.  Also ok to initialize k=0 if
  //    mode is small.
  // 2. Compute pmf(k) to high accuracy.
  // 3. Sum left-tail (<= k) likelihoods.  (Guaranteed to be < targetp unless
  //    k=0, in which case we can immediately return 0.)  Use float64 instead
  //    of dd_real precision when we can get away with it.
  // 4. Sum inward until the sum >= targetp, at which point we return k.
  const dd_real target_lnp_ddr = log_target? targetp_or_lnp_ddr : ddr_log(targetp_or_lnp_ddr);
  // Find the x on the left side where the quadratic crosses
  // y=log(targetp/mode).  If there's no such point, just start at
  // x=mode-1.
  const double target_lnprob = target_lnp_ddr.x[0] - log(mode);
  // (-b - sqrt(b^2 - 4ac)) / 2a
  const double discrim = x1_coeff * x1_coeff - 4 * x2_coeff * (x0_coeff - target_lnprob);
  double k;
  if (discrim < 0.0) {
    k = mode - 1;
  } else {
    double sqrt_discrim = sqrt(discrim);
    if (x2_coeff > 0.0) {
      // this shouldn't happen
      sqrt_discrim = -sqrt_discrim;
    }
    k = trunc(mode + (sqrt_discrim - x1_coeff) / (2 * x2_coeff));
    if (k < 0) {
      k = 0;
    }
  }
  // Use Pbinom benchmark result for now, could tune this separately later.
  const uint32_t use_bfrac = (n > 131072) && (k >= 2048);

  // |log(pmf(0) / pmf(1))| is larger than all the other gaps between adjacent
  // log(pmf()) points to the left of the mode.  So this ensures there's at
  // least one value of k where log(pmf(k)) - target_lnprob is in
  // (lnprob_diff_min, 0], letting us exit the loop; unless log(pmf(0)) >=
  // target_lnprob, in which case we exit the loop at k=0.
  double lnprob_diff_min = log(qdp_ddr.x[0] / S_CAST(double, n)) * (1 + kSmallEpsilon);
  if (!use_bfrac) {
    // Our relative error budget is usually 2^{-54}.
    // We want to ensure that accumulated error when evaluating the outer part
    // of the tailsum using plain float64 arithmetic < 2^{-55}.  Then the other
    // half of the budget covers dd_real-precision evaluation of log(pmf(k))
    // and the inner part of the tailsum.
    // 2^{-55} corresponds to 0.125 ULPs; Pbinom() comments elaborate on the
    // squared term in the denominator.  When we're using the simple tailsum
    // algorithm, there's no speedup from landing closer than tailsum_ddr_end,
    // so we widen the landing window accordingly.
    const double tailsum_ddr_end = -log(8 * mode * mode);
    if (lnprob_diff_min > tailsum_ddr_end) {
      lnprob_diff_min = tailsum_ddr_end;
    }
  }
  double nmk;
  dd_real cur_lnprob_ddr;
  while (1) {
    // probable todo: in use_bfrac case, raise target_lnprob as much as
    // possible during this loop
    nmk = n - k;
    // Evaluate pmf(k) to high precision.
    double lnprob_diff;
    if (!use_tdr) {
      cur_lnprob_ddr = ddr_sub(ddr_sort_and_add3(ddr_muld(ddr_make_td(logp_tdr), k), ddr_muld(ddr_make_td(logq_tdr), nmk), ddr_make_td(n_lfact_tdr)),
                               ddr_add_lfacts(k, nmk));
    } else {
      cur_lnprob_ddr = ddr_make_td(tdr_sub(tdr_sort_and_add3(tdr_muld(logp_tdr, k), tdr_muld(logq_tdr, nmk), n_lfact_tdr),
                                           tdr_add_lfacts(k, nmk)));
    }
    lnprob_diff = cur_lnprob_ddr.x[0] - target_lnprob;
    if (lnprob_diff > 0) {
      if (k == 0) {
        break;
      }
      // This calculation can be lower-precision.
      const double ll_deriv = log(pdq_ddr.x[0] * (nmk + 1) / k);
      k -= ceil(lnprob_diff / ll_deriv);
      if (k < 0) {
        k = 0;
      }
    } else if (lnprob_diff > lnprob_diff_min) {
      break;
    } else {
      const double ll_deriv = log(pdq_ddr.x[0] * nmk / (k + 1));
      k += S_CAST(int64_t, -lnprob_diff / ll_deriv);
    }
  }
  // Express current likelihood as a fraction of targetp.
  const double tailenter_k = k;
  const dd_real tailenter_lik_ddr = ddr_exp(ddr_sub(cur_lnprob_ddr, target_lnp_ddr));
  dd_real tailsum_ddr;
  if (use_bfrac) {
    tailsum_ddr = ddr_mul(tailenter_lik_ddr, binom_ltail_lik_bfrac_ddr(S_CAST(int64_t, k), n, succp_ddr, failp_ddr));
  } else {
    tailsum_ddr = binom_ltail_lik_simple_ddr(k, nmk, tailenter_lik_ddr, qdp_ddr, 0.125);
  }
  dd_real lik_ddr = tailenter_lik_ddr;
  while (ddr_ltd(tailsum_ddr, 1.0)) {
    k += 1;
    lik_ddr = ddr_mul(lik_ddr, ddr_divd(ddr_muld(pdq_ddr, nmk), k));
    nmk -= 1;
    tailsum_ddr = ddr_add(tailsum_ddr, lik_ddr);
  }
  return inv? (n - S_CAST(int64_t, k)) : S_CAST(int64_t, k);
}

void materialize_oddsratio_p_q_tdr(uint32_t succ_flipped, td_real* p_tdr_ptr, td_real* q_tdr_ptr, td_real* succ_odds_ratio_tdr_ptr) {
  // Currently safe to assume that either succ_odds_ratio_tdr is fully
  // evaluated, or both succ_odds_ratio_tdr and one of {p_tdr, q_tdr} needs to
  // be.
  if (succ_odds_ratio_tdr_ptr->x[2] == DBL_MAX) {
    td_real* numer_tdr_ptr = p_tdr_ptr;
    td_real* denom_tdr_ptr = q_tdr_ptr;
    if (succ_flipped) {
      numer_tdr_ptr = q_tdr_ptr;
      denom_tdr_ptr = p_tdr_ptr;
    }
    *denom_tdr_ptr = tdr_addd(tdr_negate(*numer_tdr_ptr), 1.0);
    *succ_odds_ratio_tdr_ptr = tdr_accurate_div(*numer_tdr_ptr, *denom_tdr_ptr);
  }
}

// Assumes 0 <= obs_succ <= obs_tot < 2^52 and 2^{-960} < p < 1.
//
// Note that this can be written in terms of dbinom() and pbinom():
// 1. Calculate mode, determine which tail obs_succ is on.  Early-exit if we're
//    at a mode.
// 2. Calculate pbinom() for tail, and log-dbinom() for starting table.
// 3. Search for furthest-inward value of succ on other tail where log-dbinom()
//    is <= starting log-dbinom().
// 4. Add its tail pbinom() value, return.
// Several other 2-sided exact tests based on unimodal distributions (e.g.
// Fisher's 2x2) can be performed in the same manner.
double BinomTwoSidedP(int64_t obs_succ, int64_t obs_tot, td_real p_tdr, int32_t midp, uint32_t logp) {
  if (!obs_tot) {
    if (midp) {
      return logp? -kLn2 : 0.5;
    }
    return logp? 0.0 : 1.0;
  }
  // Normalize: succ <= mode.
  // (even if there is rounding error, this is enough to guarantee that succ-1
  // has lower likelihood than succ.)
  uint32_t succ_flipped = 0;
  if (obs_succ > obs_tot * p_tdr.x[0]) {
    obs_succ = obs_tot - obs_succ;
    succ_flipped = 1;
  }
  double succ = obs_succ;
  double fail = obs_tot - obs_succ;
  // When p != 0.5, succ_odds_ratio_tdr and one of {p_tdr, q_tdr} will usually
  // only be evaluated with dd_real arithmetic (since td_real arithmetic is
  // ~5-6x as expensive).  [2] is set to DBL_MAX in this case, and
  // materialize_oddsratio_p_q_tdr() is called later as necessary.
  const uint32_t p_is_half = tdr_is(p_tdr, 0.5);
  td_real q_tdr;
  td_real succ_odds_ratio_tdr;
  double first_inward_mult;
  if (p_is_half) {
    if ((!midp) && (fail <= succ + 1)) {
      // could use binom_ln_prob_approx() to accelerate midp subcase
      return logp? 0.0 : 1.0;
    }
    q_tdr = p_tdr;
    succ_odds_ratio_tdr = tdr_make1(1);
    first_inward_mult = fail / (succ + 1);
  } else {
    dd_real p_ddr = ddr_make_td(p_tdr);
    // p_tdr in [1 - 2^{-54}, 1) is a potentially annoying edge case.
    // Fortunately, even when the difference from 1 is much smaller than
    // 2^{-54}, p_ddr.x[1] still accurately represents the difference from 1.
    dd_real q_ddr = ddr_negate(ddr_make_td(tdr_addd(p_tdr, -1.0)));
    if (succ_flipped) {
      swap_ddr(&p_ddr, &q_ddr);
    }
    dd_real succ_odds_ratio_ddr = ddr_accurate_div(p_ddr, q_ddr);
    first_inward_mult = fail * succ_odds_ratio_ddr.x[0] / (succ + 1);
    if ((!midp) && (first_inward_mult <= 1 + 2 * k2m52)) {
      if (first_inward_mult < 1 - 2 * k2m52) {
        return logp? 0.0 : 1.0;
      }
      q_tdr = tdr_addd(tdr_negate(p_tdr), 1.0);
      if (succ_flipped) {
        swap_tdr(&p_tdr, &q_tdr);
      }
      succ_odds_ratio_tdr = tdr_accurate_div(p_tdr, q_tdr);
      const td_real first_inward_mult_numer_tdr = tdr_muld(succ_odds_ratio_tdr, fail);
      if (tdr_leq(first_inward_mult_numer_tdr, tdr_make2(succ + 1, (succ + 1) * (4 * _tdr_eps)))) {
        return logp? 0.0 : 1.0;
      }
      // no need to update first_inward_mult
    } else {
      succ_odds_ratio_tdr = tdr_make(succ_odds_ratio_ddr.x[0], succ_odds_ratio_ddr.x[1], DBL_MAX);
      if (!succ_flipped) {
        q_tdr = tdr_make(q_ddr.x[0], q_ddr.x[1], DBL_MAX);
      } else {
        q_tdr = p_tdr;
        p_tdr = tdr_make(p_ddr.x[0], p_ddr.x[1], DBL_MAX);
      }
    }
  }
  const double succ_odds_ratio = succ_odds_ratio_tdr.x[0];
  // todo: benchmark different thresholds
  // const uint32_t use_bfrac = (obs_tot > 512) && (MINV(obs_succ + 1, obs_tot - obs_succ) >= 40);
  double tail_sum = binom_ltail_lik_simple(succ, fail, succ_odds_ratio, midp);
  /*
  double tail_sum;
  if (use_bfrac) {
    tail_sum = binom_tail_lik_bfrac(obs_succ, obs_tot, ddr_make_td(p_tdr), ddr_make_td(q_tdr), 0, midp);
  } else {
    tail_sum = binom_ltail_lik_simple(succ, fail, succ_odds_ratio, midp);
  }
  */
  // In the common case, where we're close enough to the mode that float64
  // underflow/overflow isn't an issue, use the original algorithm: sum all
  // center relative-likelihoods, sum far-tail relative-likelihoods to
  // floating-point precision limit, return
  //   tail_sum / (tail_sum + center_sum)
  // or its log.
  //
  // Unfortunately, an extremal rate makes center_sum overflow possible much
  // closer to the mode than is the case for the Fisher/HWE exact tests.  So
  // instead of checking whether we're a constant number of steps from the
  // mode, we compute a lower bound on the number of non-overflowing inward
  // steps we can take, using the log of the first-inward-step multiplier
  // (subsequent steps have smaller multipliers).
  succ = obs_succ;
  fail = obs_tot - obs_succ;
  const double ln_mult = log(first_inward_mult);
  double overflow_steps_lower_bound = 0x7fffffff;
  // log(DBL_MAX / (2^31 - 1)) = 688.295...
  // probably want a different rule in use_bfrac case
  if (ln_mult > (688.295 / S_CAST(double, 0x7fffffff))) {
    overflow_steps_lower_bound = 688.295 / ln_mult;
  }
  // possible for modal_succ to round up to just obs_tot
  const double obs_totd = obs_tot;
  const double modal_succ = obs_totd * p_tdr.x[0];
  if (succ + overflow_steps_lower_bound > modal_succ) {
    double one_plus_scaled_eps = 1 + k2m52;
    double center_sum = midp * 0.5;
    double lik = 1;
    while (1) {
      succ += 1;
      lik *= succ_odds_ratio * fail / succ;
      fail -= 1;
      // succ_odds_ratio is off by up to 0.5 ULP, and we have two multiplies
      // and a divide.
      one_plus_scaled_eps += 2 * k2m52;
      if (lik < one_plus_scaled_eps) {
        if (lik <= 2 - one_plus_scaled_eps) {
          tail_sum += lik;
          break;
        }
        // Near-tie.  True value of lik can be greater than, equal to, or
        // less than 1.
        materialize_oddsratio_p_q_tdr(succ_flipped, &p_tdr, &q_tdr, &succ_odds_ratio_tdr);
        td_real starting_lnprobv_tdr = tdr_make1(DBL_MAX);
        td_real ln_odds_ratio_tdr = tdr_make1(DBL_MAX);
        const intptr_t cmp_result = BinomCompare(obs_succ, obs_tot, succ_odds_ratio_tdr, S_CAST(int64_t, succ), &starting_lnprobv_tdr, &ln_odds_ratio_tdr, &lik);
        one_plus_scaled_eps = 1 + 3 * k2m52;
        if (cmp_result <= 0) {
          tail_sum += lik;
          if (midp && (cmp_result == 0)) {
            tail_sum -= 0.5;
            center_sum += 0.5;
          }
          break;
        }
      }
      center_sum += lik;
    }
    // Continue down tail to floating-point precision limit.
    while (1) {
      succ += 1;
      lik *= succ_odds_ratio * fail / succ;
      fail -= 1;
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
    }
    const double pval = tail_sum / (tail_sum + center_sum);
    return logp? log(pval) : pval;
  }
  const uint32_t tdr_lnprobv_needed = use_tdr_for_binom_lnprob(obs_tot);
  td_real ln_odds_ratio_tdr;
  td_real starting_lnprobv_tdr;
  dd_real starting_lnprob_ddr;
  if (!tdr_lnprobv_needed) {
    const dd_real ln_odds_ratio_ddr = ddr_log(ddr_make_td(succ_odds_ratio_tdr));
    const dd_real starting_lnprobv_ddr =
      ddr_sub(ddr_muld(ln_odds_ratio_ddr, succ),
              ddr_add_lfacts(succ, fail));
    dd_real lnfail_ddr = _ddr_log05;
    if (!p_is_half) {
      if (!succ_flipped) {
        lnfail_ddr = ddr_log1p(ddr_negate(ddr_make_td(p_tdr)));
      } else {
        lnfail_ddr = ddr_log(ddr_make_td(q_tdr));
      }
    }
    const dd_real lnprobf_ddr =
      ddr_add(ddr_lfact(obs_totd), ddr_muld(lnfail_ddr, obs_totd));
    starting_lnprob_ddr = ddr_add(lnprobf_ddr, starting_lnprobv_ddr);
    starting_lnprobv_tdr = tdr_make(starting_lnprobv_ddr.x[0], starting_lnprobv_ddr.x[1], DBL_MAX);
    if (p_is_half) {
      ln_odds_ratio_tdr = tdr_make1(0);
    } else {
      ln_odds_ratio_tdr = tdr_make(ln_odds_ratio_ddr.x[0], ln_odds_ratio_ddr.x[1], DBL_MAX);
    }
  } else {
    materialize_oddsratio_p_q_tdr(succ_flipped, &succ_odds_ratio_tdr, &p_tdr, &q_tdr);
    ln_odds_ratio_tdr = tdr_log(succ_odds_ratio_tdr);
    starting_lnprobv_tdr =
      tdr_sub(tdr_muld(ln_odds_ratio_tdr, succ),
              tdr_add_lfacts(succ, fail));
    td_real lnfail_tdr = _tdr_log05;
    if (!p_is_half) {
      // probable todo: this is a bit redundant with earlier initialization
      if (!succ_flipped) {
        // handle tiny p correctly
        lnfail_tdr = tdr_log1p(tdr_negate(p_tdr));
      } else {
        lnfail_tdr = tdr_log(q_tdr);
      }
    }
    const td_real lnprobf_tdr =
      tdr_add(tdr_lfact(obs_totd), tdr_muld(lnfail_tdr, obs_totd));
    starting_lnprob_ddr = ddr_make_td(tdr_add(lnprobf_tdr, starting_lnprobv_tdr));
  }

  // Now we want to jump near the other tail, without evaluating that many
  // contingency table log-likelihoods along the way.
  //
  // Each full log-likelihood evaluation requires 2 ddr_lfact() or tdr_lfact()
  // calls.  Since they are now performed with extra precision, they require
  // hundreds or thousands of floating-point operations, so we want to limit
  // ourselves to 1-2 full evaluations most of the time.  (Possible todo: use
  // lower-accuracy Lfact() to jump around, followed by {d,q}dr_lfact() when
  // exiting the loop.  Should be an easy performance win, but there's a
  // complexity cost so I'll wait until I see a scenario where this branch
  // executes frequently...)
  //
  // The current heuristic starts by reflecting (succ - 1) across the
  // (continuous) mode, performing a full log-likelihood check at the nearest
  // valid point.  Hopefully we find that we're in (starting_lnprob -
  // tolerance, starting_lnprob], so we're at or near a table that contributes
  // non-negligibly to the tail-sum; unlike the Fisher's and HWE cases, we
  // can't fix the tolerance at 62 * kLn2, but we can compute a value >= 53 *
  // kLn2 large enough to guarantee at least one point falls inside.
  //
  // If not, we jump again, using Newton's method.
  // If succ is too low (i.e. current log-likelihood is too high), when we
  // increase succ by 1, the likelihood gets multiplied by
  //   succ_odds_ratio * fail / (succ+1)
  // i.e. we're adding the logarithm of this value to the log-likelihood.
  // If succ is too high, when we decrease succ by 1, the likelihood gets
  // multiplied by
  //   succ / (succ_odds_ratio * (fail+1))
  // We use the log of the first expression as the Newton's method f'(x) when
  // we're jumping to higher succ, and the negative-log of the second
  // expression when we're jumping to lower homr.
  // f''(x) is always negative, so we can aim for starting_lnprob instead of
  // the middle of the interval.

  // L(obs_tot) / L(obs_tot-1) = succ_odds_ratio * 1 / obs_tot
  // If this value is >= 1, obs_tot is a mode.  Separating out that case makes
  // the remaining logic simpler.
  if (tdr_geq(succ_odds_ratio_tdr, tdr_make2(obs_totd, obs_totd * (4 * _tdr_eps)))) {
    return join_log_and_nonlog(starting_lnprob_ddr, tail_sum, logp);
  }

  succ = 2 * modal_succ - (succ - 1);
  if (succ > obs_totd) {
    succ = obs_totd;
  }
  succ = S_CAST(int64_t, succ);

  // obs_tot is past the mode, and |log(L(obs_tot) / L(obs_tot-1))| is the
  // largest gap between adjacent log-likelihoods on this tail.  Set
  // |lnprobv_diff_min| >= this value.
  double lnprobv_diff_min = log(succ_odds_ratio / obs_totd) * (1 + kSmallEpsilon);
  // if (use_bfrac && (lnprobv_diff_min > -53 * kLn2)) {
  if (lnprobv_diff_min > -53 * kLn2) {
    lnprobv_diff_min = -53 * kLn2;
  }

  double lik;
  while (1) {
    fail = obs_totd - succ;
    double lnprobv_diff;
    if (!tdr_lnprobv_needed) {
      const dd_real lnprobv_ddr =
        ddr_sub(ddr_muld(ddr_make_td(ln_odds_ratio_tdr), succ),
                ddr_add_lfacts(succ, fail));
      lnprobv_diff = ddr_sub(lnprobv_ddr, ddr_make_td(starting_lnprobv_tdr)).x[0];
    } else {
      const td_real lnprobv_tdr =
        tdr_sub(tdr_muld(ln_odds_ratio_tdr, succ),
                tdr_add_lfacts(succ, fail));
      lnprobv_diff = tdr_sub(lnprobv_tdr, starting_lnprobv_tdr).x[0];
    }
    if (lnprobv_diff >= k2m53) {
      if (fail == 0) {
        return join_log_and_nonlog(starting_lnprob_ddr, tail_sum, logp);
      }
      const double ll_deriv = ln_odds_ratio_tdr.x[0] + log(fail / (succ + 1));
      // Round up, to guarantee that we make progress.
      // (lnprobv_diff is positive and ll_deriv is negative.)
      // This may overshoot.  But the function is guaranteed to terminate
      // because we never overshoot (and we do always make progress on each
      // step) once we're on the other side.
      succ += ceil(-lnprobv_diff / ll_deriv);
      if (succ > obs_totd) {
        succ = obs_totd;
      }
    } else if (lnprobv_diff > lnprobv_diff_min) {
      lik = exp(lnprobv_diff);
      break;
    } else {
      const double ll_deriv = ln_odds_ratio_tdr.x[0] + log((fail + 1) / succ);
      // Round down, to guarantee we don't overshoot.
      // |lnprobv_diff| >= |lnprobv_diff_min| > |ll_deriv| so we're guaranteed
      // to make progress.
      succ -= S_CAST(int64_t, lnprobv_diff / ll_deriv);
      // todo: if use_bfrac, tighten lnprobv_diff_min since we are searching
      // for the exact crossover point.
    }
  }
  // Sum toward center, until lik >= 1.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  // Save where we're starting on this tail, which isn't necessarily on the
  // boundary.  We sum inward until relative-likelihood > 1, then we jump back
  // to tailenter_succ and sum outward.
  const double tailenter_lik = lik;
  const double tailenter_succ = succ;
  while (lik <= one_minus_scaled_eps) {
    tail_sum += lik;
    fail += 1;
    lik *= succ / (succ_odds_ratio * fail);
    succ -= 1;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lik < 2 - one_minus_scaled_eps) {
    materialize_oddsratio_p_q_tdr(succ_flipped, &succ_odds_ratio_tdr, &p_tdr, &q_tdr);
    const intptr_t cmp_result = BinomCompare(obs_succ, obs_tot, succ_odds_ratio_tdr, S_CAST(int64_t, succ), &starting_lnprobv_tdr, &ln_odds_ratio_tdr, &lik);
    if (cmp_result <= 0) {
      tail_sum += lik;
      if (midp && (cmp_result == 0)) {
        tail_sum -= 0.5;
      }
    }
  }
  // Sum away from center, until sums stop changing.
  lik = tailenter_lik;
  succ = tailenter_succ;
  fail = obs_totd - succ;
  while (1) {
    succ += 1;
    lik *= succ_odds_ratio * fail / succ;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
    fail -= 1;
  }
  return join_log_and_nonlog(starting_lnprob_ddr, tail_sum, logp);
}


#ifdef __cplusplus
}
#endif
