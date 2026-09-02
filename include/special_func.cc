// Copyright (C) 2026 Christopher Chang.
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

#include "special_func.h"

#include <assert.h>
#include <math.h>

#include "plink2_base.h"
#include "plink2_float.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// float64- (and a few higher-precision) routines for evaluating special
// functions.  Currently touches beta and inverse-error.

// binom_ln_prob_loader() implements Catherine Loader's algorithm:
//   https://www.r-project.org/doc/reports/CLoader-dbinom-2002.pdf
// with dd_reals (similar in character to R ebd0()).  Relative error should be
// better than ~2^{-90}?
//
// The key idea is to decompose the log-probability into a nonpositive term
// corresponding to p_0 := n/k, and two more nonpositive terms of the form
//   C * (x log x + 1 - x).
// Since (x log x + 1 - x) can be accurately evaluated via series expansion for
// x near 1, we never have significant cancellation.
//
// Some of the overflow-avoidance logic is derived from GPL-2 R code.
dd_real loader_bd0(dd_real x_ddr, dd_real np_ddr) {
  const dd_real x_minus_np_ddr = ddr_sub(x_ddr, np_ddr);
  // x+np may overflow.
  const dd_real half_x_plus_np_ddr = ddr_add(ddr_mul_pwr2(x_ddr, 0.5), ddr_mul_pwr2(np_ddr, 0.5));
  if (fabs(x_minus_np_ddr.x[0]) >= 0.2 * half_x_plus_np_ddr.x[0]) {
    dd_real log_x_div_np_ddr;
    // Avoid potential x/np overflow.
    if (x_ddr.x[0] * (1.0 / (k2p800 * k2p100)) < np_ddr.x[0]) {
      log_x_div_np_ddr = ddr_log(ddr_accurate_div(x_ddr, np_ddr));
    } else {
      log_x_div_np_ddr = ddr_sub(ddr_log(x_ddr), ddr_log_extdomain(np_ddr));
    }
    // Avoid potential ddr_mul() overflow.
    return ddr_mul_pwr2(ddr_sub(ddr_mul(ddr_mul_pwr2(x_ddr, 1.0 / 2048), log_x_div_np_ddr), ddr_mul_pwr2(x_minus_np_ddr, 1.0 / 2048)), 2048);
  }
  const dd_real double_v_ddr = ddr_accurate_div(x_minus_np_ddr, half_x_plus_np_ddr);
  dd_real ej_ddr = ddr_mul(x_ddr, double_v_ddr);
  const dd_real v_ddr = ddr_mul_pwr2(double_v_ddr, 0.5);
  dd_real s_ddr = ddr_mul(x_minus_np_ddr, v_ddr);
  const dd_real v2_ddr = ddr_sqr(v_ddr);
  for (double j = 1; ; j += 1) {
    ej_ddr = ddr_mul(ej_ddr, v2_ddr);
    const dd_real s1_ddr = ddr_add(s_ddr, ddr_divd(ej_ddr, 2*j + 1));
    if (ddr_eq(s1_ddr, s_ddr)) {
      return s_ddr;
    }
    s_ddr = s1_ddr;
  }
}

void binom_ln_prob_loader_part1(dd_real k_ddr, dd_real nmk_ddr, dd_real n_ddr, dd_real* stirlerr_ddr_ptr, dd_real* half_lf_ddr_ptr) {
  const dd_real stirlerr_1_ddr = ddr_stirlerr(n_ddr);
  dd_real stirlerr_2_ddr = ddr_stirlerr(k_ddr);
  dd_real stirlerr_3_ddr = ddr_stirlerr(nmk_ddr);
  if (stirlerr_2_ddr.x[0] > stirlerr_3_ddr.x[0]) {
    swap_ddr(&stirlerr_2_ddr, &stirlerr_3_ddr);
  }
  *stirlerr_ddr_ptr = ddr_sub(ddr_sub(stirlerr_1_ddr, stirlerr_2_ddr), stirlerr_3_ddr);
  // Avoid potential overflow/underflow in Loader's original code.  See R
  // src/nmath/dbinom.c .
  const dd_real lf_ddr =
    ddr_add3(ddr_log1p(ddr_accurate_div(ddr_negate(k_ddr),
                                        n_ddr)),
             ddr_mul_pwr2(_ddr_half_log_2pi, 2),
             ddr_log(k_ddr));
  *half_lf_ddr_ptr = ddr_mul_pwr2(lf_ddr, 0.5);
}

dd_real binom_ln_prob_loader_part2(dd_real k_ddr, dd_real nmk_ddr, dd_real n_ddr, dd_real p_ddr, dd_real q_ddr, dd_real stirlerr_ddr, dd_real half_lf_ddr) {
  dd_real ddrs[3];
  ddrs[0] = stirlerr_ddr;
  ddrs[1] = ddr_negate(loader_bd0(k_ddr, ddr_mul(n_ddr, p_ddr)));
  ddrs[2] = ddr_negate(loader_bd0(nmk_ddr, ddr_mul(n_ddr, q_ddr)));
  const dd_real lc_ddr = ddr_sort_and_add(3, ddrs);
  return ddr_sub(lc_ddr, half_lf_ddr);
}

dd_real binom_ln_prob_loader(dd_real k_ddr, dd_real n_ddr, dd_real p_ddr, dd_real q_ddr) {
  // Assumes k <= n are nonnegative integers where ddr_sub(n_ddr, k_ddr)
  // does not incur any error in representing n-k.
  // Assumes 0 < p,q < 1, p+q=1; one of them may be denormal.
  if (ddr_is_zero(k_ddr)) {
    return ddr_mul(ddr_log_extdomain_maybehalf(q_ddr), n_ddr);
  }
  const dd_real nmk_ddr = ddr_sub(n_ddr, k_ddr);
  if (ddr_is_zero(nmk_ddr)) {
    return ddr_mul(ddr_log_extdomain_maybehalf(p_ddr), n_ddr);
  }
  dd_real stirlerr_ddr;
  dd_real half_lf_ddr;
  binom_ln_prob_loader_part1(k_ddr, nmk_ddr, n_ddr, &stirlerr_ddr, &half_lf_ddr);
  return binom_ln_prob_loader_part2(k_ddr, nmk_ddr, n_ddr, p_ddr, q_ddr, stirlerr_ddr, half_lf_ddr);
}

// ibeta_...() and dependencies below are adapted from Boost 1.91.0.
// This derived code is subject to the following license:
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

// this depends on the polynomial coefficients above
// exactly 808618867 * 2^{-27}, don't need to represent this as dd_real
static const double kLanczosDoubleG = 6.024680040776729583740234375;

double lanczos_sum_d_expg_scaled_imp(double zz, double* s2_ptr) {
  // zz currently guaranteed to be >1.
  zz = 1 / zz;
  const double s1 = POLY12(zz,
                           0.006061842346248907,
                           0.5098416655656676,
                           19.519927882476175,
                           449.9445569063168,
                           6955.999602515376,
                           75999.29304014542,
                           601859.6171681099,
                           3481712.154980646,
                           14605578.087685067,
                           43338889.32467614,
                           86363131.2881386,
                           103794043.11634454,
                           56906521.913471565);
  *s2_ptr = POLY11(zz,
                   1,
                   66,
                   1925,
                   32670,
                   357423,
                   2637558,
                   13339535,
                   45995730,
                   105258076,
                   150917976,
                   120543840,
                   39916800);
  return s1;
}

dd_real ibeta_power_terms_d_ln(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr) {
  // Returns log((p^a)(q^b) / Beta(a,b))
  //       = log((p^a)(q^b)(a+b-1)! / ((a-1)!(b-1)!)).
  //
  // todo: compare to binom_ln_prob_loader()-based approach
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
  const dd_real term1_ddr =
    ddr_sqr(
      ddr_accurate_div(ddr_muld(ddr_mul2d(numer_a, numer_b), numer_c),
                       ddr_muld(ddr_mul2d(denom_a, denom_b), denom_c)));
  // now multiplies by bgh later to avoid potential overflow
  const dd_real term2_ddr =
    ddr_mul(ddr_accurate_div(agh_ddr,
                             ddr_mul(cgh_ddr, _ddr_e)),
            bgh_ddr);
  dd_real result_ddr =
    ddr_mul_pwr2(ddr_log(ddr_mul(term1_ddr, term2_ddr)),
                 0.5);

  // Calculate l1 and l2 with extra precision, since magnitude can greatly
  // exceed that of ln(result).
  // This removes the need for special cases.
  const dd_real l1_ddr =
    ddr_accurate_div(ddr_negate(ddr_add(aq_minus_bp_ddr,
                                        ddr_muld(q_ddr, gh))),
                     agh_ddr);
  const dd_real l2_ddr =
    ddr_accurate_div(ddr_sub(aq_minus_bp_ddr,
                             ddr_muld(p_ddr, gh)),
                     bgh_ddr);
  return ddr_sort_and_add3(result_ddr,
                           ddr_muld(ddr_log1p(l1_ddr), aa),
                           ddr_muld(ddr_log1p(l2_ddr), bb));
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
    const double shared_frac = (bb - mm) / denom;
    const double cur_a = ((aa + mm - 1) / denom) * ((aa + bb + mm - 1) * xx * shared_frac) * (mm * xx);
    const double cur_b = prefer_fma((aa + mm) / (aa + 2 * mm + 1), prefer_fma(mm, two_minus_x, ay_minus_bx_plus1), mm + shared_frac * mm * xx);
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

/*
double erfcx(double x) {
  // e^(x^2) * erfc(x).
  //
  // This is based on Norbert Juffa's implementation from
  // https://stackoverflow.com/a/39777361 , which adapts the approach in
  //   Shepherd MM and Laframboise JG (1981) Chebyshev Approximation of (1+2x)
  //   exp(x^2) erfc x in 0 <= x < \inf .  Mathematics of Computation, 36.
  // to achieve ~3 ULP maximum error with float64 arithmetic.  CC BY-SA 3.0
  // license.
  //
  // This function makes heavy use of FMA to limit rounding error, so it
  // incurs a significant performance hit when compiled to support older x86.
  // (We could use prefer_fma() in place of fma(), but Juffa's maximum-error
  // measurement was with FMA so I'd rather not risk it.)
  //
  // An alternative is the erfcx_y100() function from
  // http://ab-initio.mit.edu/faddeeva/ , which uses a huge switch statement
  // with a different degree-6 polynomial for each of 100 cases.  That might
  // make sense for an application where erfcx() was the primary compute
  // bottleneck, but it's only secondary for BASYM so I doubt it's worth the
  // i-cache hit.

  const double a = fabs(x);  // no need to preserve NaN type

  // Compute q = (a-4)/(a+4) accurately.  [0, \inf) -> [-1, 1].
  // (This looks awkward, but compare to
  //   ddr_accurate_div(ddr_add2d(a, -4.0), ddr_add2d(a, 4.0)).x[0]
  // ...
  // todo: understand why this works.)
  double q;
  {
    const double am4 = a - 4.0;
    const double ap4 = a + 4.0;
    const double recip_ap4 = 1.0 / ap4;
    const double q_approx = am4 * recip_ap4;
    const double negt = 4 * (q_approx + 1.0) - a;
    const double e = fma(q_approx, -a, -negt);
    q = fma(recip_ap4, e, q_approx);
  }

  // Approximate (1+2*a)*exp(a*a)*erfc(a) as p(q)+1 for q in [-1, 1].
  // This currently deviates slightly from Juffa in using 4 FMA chains instead
  // of 2; could change that if it matters.
  const double p = POLY23_FMA_SOLO(q,
                                   2.3299511862555250e-01,
                                   -1.3962111684056208e-01,
                                   1.5379652102610957e-02,
                                   6.8097054254651804e-02,
                                   -1.0103906603588378e-01,
                                   9.3732834999538536e-02,
                                   -6.6330365820039094e-02,
                                   3.7167515521269866e-02,
                                   -1.6197733983519948e-02,
                                   5.0319701025945277e-03,
                                   -7.5777369791018515e-04,
                                   -1.9925728768782324e-04,
                                   1.5062307184282616e-04,
                                   -2.4397380523258482e-05,
                                   -1.1225056665965572e-05,
                                   5.7059822144459833e-06,
                                   2.9796165315625938e-07,
                                   -8.2040389712752056e-07,
                                   7.1190423171700940e-08,
                                   1.0585794011876720e-07,
                                   -1.6386753783877791e-08,
                                   -1.2155985739342269e-08,
                                   1.5764464777959401e-09,
                                   8.9820305531190140e-10);

  // Divide (1+p) by (1+2*a) -> exp(a*a)*erfc(a)
  double result;
  {
    const double d = a + 0.5;
    // Juffa separated out a (1.0 / d) reciprocal operation, but I don't think
    // that's relevant to the target x86/ARM platforms?
    const double half_recip_d = 0.5 / d;
    q = fma(p, half_recip_d, half_recip_d);  // q = (p+1)/(1+2*a)
    const double t = 2 * q;
    const double e = (p - q) + fma(t, -a, 1.0);  // residual: (p+1)-q*(1+2*a)
    result = fma(e, half_recip_d, q);
  }
  // Handle negative arguments: erfcx(x) = 2*exp(x*x) - erfcx(|x|)
  if (x < 0.0) {
    const double s = x * x;
    const double d = fma(x, x, -s);
    const double e = exp(s);
    result = e - result;
    result = fma(e, 2 * d, result);
    result = result + e;
    if (e > DBL_MAX) {
      result = e;  // avoid creating NaN
    }
  }
  return result;
}
*/

dd_real erfcx_ddr(dd_real x_ddr) {
  // e^(x^2) * erfc(x).

  // We want a few more bits of accuracy than TOMS 708 erfc1() provides, and
  // don't need bleeding-edge speed for now (just need to be better than
  // evaluating >100000 continued-fraction terms...).  Hopefully this port of
  // the Netlib CALERF function, which specifies constants to 18 digits and is
  // supposed to have maximum relative error slightly better than 10^{-18}, is
  // good enough.
  // If not:
  // - The Shepherd/Laframboise paper appears to have constants which yield
  //   ~22-digit accuracy.
  // - Zaghloul MR (2022) "Efficient multiple-precision computation of the
  //   scaled complementary error function and the Dawson Integral" looks like
  //   a description of an implementation with ~30-digit accuracy, but I'm not
  //   aware of a permissively-licensed implementation.

  const dd_real y_ddr = (x_ddr.x[0] < 0)? ddr_negate(x_ddr) : x_ddr;
  if (y_ddr.x[0] <= 0.46875) {
    const dd_real ysq_ddr = ddr_sqr(x_ddr);
    const dd_real a_ddr[5] = {
      {3.1611237438705655, 1.0548186774540226e-16},
      {113.86415415105016, -2.2479264470748604e-15},
      {377.485237685302, -3.6016967575997115e-15},
      {3209.3775891384694, 9.497032806277276e-14},
      {0.18577770618460315, 5.416702813818119e-21}
    };
    const dd_real b_ddr[4] = {
      {23.601290952344122, -1.431396647496149e-15},
      {244.02463793444417, 9.395336285233499e-16},
      {1282.6165260773723, -4.756335288286209e-14},
      {2844.236833439171, -2.0311236299574375e-13}
    };
    dd_real xnum_ddr = ddr_mul(ysq_ddr, a_ddr[4]);
    dd_real xden_ddr = ysq_ddr;
    for (uint32_t i = 0; i < 3; ++i) {
      xnum_ddr = ddr_mul(ddr_add(xnum_ddr, a_ddr[i]), ysq_ddr);
      xden_ddr = ddr_mul(ddr_add(xden_ddr, b_ddr[i]), ysq_ddr);
    }
    const dd_real result_ddr = ddr_accurate_div(ddr_mul(x_ddr, ddr_add(xnum_ddr, a_ddr[3])), ddr_add(xden_ddr, b_ddr[3]));
    return ddr_mul(ddr_exp(ysq_ddr), ddr_negate(ddr_addd(result_ddr, -1)));
  }
  dd_real result_ddr;
  if (y_ddr.x[0] <= 4.0) {
    const dd_real c_ddr[9] = {
      {0.5641884969886701, -4.053940047736978e-17},
      {8.883149794388377, -8.786879024142399e-16},
      {66.11919063714163, 1.4691960141062737e-15},
      {298.6351381974001, 6.773012172430754e-15},
      {881.952221241769, 3.7460106108337643e-14},
      {1712.0476126340707, -9.41757458075881e-14},
      {2051.0783778260716, -1.0350867055356502e-13},
      {1230.3393547979972, 4.730488546192646e-14},
      {2.1531153547440383e-08, 1.1924348055929413e-24}
    };
    const dd_real d_ddr[8] = {
      {15.744926110709835, -7.983691024244763e-16},
      {117.6939508913125, 2.235159965697676e-15},
      {537.1811018620099, -2.409748174995184e-14},
      {1621.3895745666903, -8.124007698148489e-14},
      {3290.7992357334597, -5.1796375662088396e-14},
      {4362.619090143247, -2.1654533237218857e-13},
      {3439.3676741437216, 2.4890174493193628e-14},
      {1230.3393548037495, -1.093101528659463e-13}
    };
    dd_real xnum_ddr = ddr_mul(c_ddr[8], y_ddr);
    dd_real xden_ddr = y_ddr;
    for (uint32_t i = 0; i < 7; ++i) {
      xnum_ddr = ddr_mul(ddr_add(xnum_ddr, c_ddr[i]), y_ddr);
      xden_ddr = ddr_mul(ddr_add(xden_ddr, d_ddr[i]), y_ddr);
    }
    result_ddr = ddr_accurate_div(ddr_add(xnum_ddr, c_ddr[7]), ddr_add(xden_ddr, d_ddr[7]));
  } else {
    result_ddr = ddr_maked(0);
    const dd_real _ddr_inv_sqrt_pi = {{0.5641895835477563, 7.66772980658294e-18}};
    if (y_ddr.x[0] >= 6.71e7) {
      if (y_ddr.x[0] < 2.53e307) {
        result_ddr = ddr_accurate_div(_ddr_inv_sqrt_pi, y_ddr);
      }
    } else {
      const dd_real p_ddr[6] = {
        {0.30532663496123236, -1.3134974790824345e-17},
        {0.36034489994980445, -1.459828991058748e-17},
        {0.12578172611122926, -1.2391333019062585e-17},
        {0.016083785148742275, 1.4105070992778564e-18},
        {0.0006587491615298378, -1.9665056634718896e-20},
        {0.016315387137302097, 3.855181504311986e-19}
      };
      const dd_real q_ddr[5] = {
        {2.568520192289822, 2.1844100914488081e-16},
        {1.8729528499234604, 2.0056509142450522e-17},
        {0.5279051029514285, -3.898581574347918e-17},
        {0.06051834131244132, -4.6397944450654906e-20},
        {0.0023352049762686918, 8.153490887023054e-20}
      };
      dd_real ysq_ddr = ddr_accurate_div(ddr_maked(1), ddr_sqr(y_ddr));
      dd_real xnum_ddr = ddr_mul(p_ddr[5], ysq_ddr);
      dd_real xden_ddr = ysq_ddr;
      for (uint32_t i = 0; i < 4; ++i) {
        xnum_ddr = ddr_mul(ddr_add(xnum_ddr, p_ddr[i]), ysq_ddr);
        xden_ddr = ddr_mul(ddr_add(xden_ddr, q_ddr[i]), ysq_ddr);
      }
      result_ddr = ddr_accurate_div(ddr_mul(ysq_ddr, ddr_add(xnum_ddr, p_ddr[4])), ddr_add(xden_ddr, q_ddr[4]));
      result_ddr = ddr_accurate_div(ddr_sub(_ddr_inv_sqrt_pi, result_ddr), y_ddr);
    }
  }
  if (x_ddr.x[0] < 0) {
    // don't need this for pbinom(), but it's nice to provide an erfcx()
    // library function
    const dd_real exponential_term_ddr = ddr_mul_pwr2(ddr_exp(ddr_sqr(x_ddr)), 2);
    if (exponential_term_ddr.x[0] == INFINITY_D) {
      // ddr_sub degrades infinity to NaN
      result_ddr = exponential_term_ddr;
    } else {
      result_ddr = ddr_sub(exponential_term_ddr, result_ddr);
    }
  }
  return result_ddr;
}

// log1pmx() and supporting logcf() are from R src/nmath/pgamma.c .
// Some recent discussion at
//   https://cran.r-project.org/web/packages/DPQ/vignettes/log1pmx-etc.pdf
// ; but I don't see a separate log1pmx() implementation in DPQ 0.6-1 so I'm
// guessing the algorithm is best left alone for now (outside of dd_real
// widening).
dd_real ddr_logcf_d2(dd_real x_ddr, double i, double eps) {
  const double d = 2;
  // Continued fraction for calculation of
  //   1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ...
  // d currently assumed to be a power of 2, update b2_ddr initialization if
  // that can no longer be assumed.
  const double scalefactor = k2p200 * (1LL << 56);
  double c1 = 2 * d;
  double c2 = i + d;
  double c4 = c2 + d;
  dd_real a1_ddr = ddr_maked(c2);
  dd_real b1_ddr = ddr_muld(ddr_addd(ddr_muld(x_ddr, -i), c2), i);
  dd_real b2_ddr = ddr_mul_pwr2(x_ddr, d*d);
  dd_real a2_ddr = ddr_sub(ddr_mul2d(c4, c2), b2_ddr);

  b2_ddr = ddr_sub(ddr_muld(b1_ddr, c4), ddr_muld(b2_ddr, i));

  while (fabs(ddr_sub(ddr_mul(a2_ddr, b1_ddr), ddr_mul(a1_ddr, b2_ddr)).x[0]) >
         fabs(eps * b1_ddr.x[0] * b2_ddr.x[0])) {
    dd_real c3_ddr = ddr_muld(x_ddr, c2 * c2);
    c2 += d;
    c4 += d;
    a1_ddr = ddr_sub(ddr_muld(a2_ddr, c4), ddr_mul(a1_ddr, c3_ddr));
    b1_ddr = ddr_sub(ddr_muld(b2_ddr, c4), ddr_mul(b1_ddr, c3_ddr));

    c3_ddr = ddr_muld(x_ddr, c1 * c1);
    c1 += d;
    c4 += d;
    a2_ddr = ddr_sub(ddr_muld(a1_ddr, c4), ddr_mul(a2_ddr, c3_ddr));
    b2_ddr = ddr_sub(ddr_muld(b1_ddr, c4), ddr_mul(b2_ddr, c3_ddr));

    if (fabs(b2_ddr.x[0]) > scalefactor) {
      a1_ddr = ddr_mul_pwr2(a1_ddr, 1 / scalefactor);
      b1_ddr = ddr_mul_pwr2(b1_ddr, 1 / scalefactor);
      a2_ddr = ddr_mul_pwr2(a2_ddr, 1 / scalefactor);
      b2_ddr = ddr_mul_pwr2(b2_ddr, 1 / scalefactor);
    } else if (fabs(b2_ddr.x[0]) < 1 / scalefactor) {
      a1_ddr = ddr_mul_pwr2(a1_ddr, scalefactor);
      b1_ddr = ddr_mul_pwr2(b1_ddr, scalefactor);
      a2_ddr = ddr_mul_pwr2(a2_ddr, scalefactor);
      b2_ddr = ddr_mul_pwr2(b2_ddr, scalefactor);
    }
  }
  return ddr_accurate_div(a2_ddr, b2_ddr);
}

static const dd_real _ddr_2_3rds = {{6.6666666666666662966e-01, 3.7007434154171882626e-17}};
static const dd_real _ddr_sqrt2 = {{1.4142135623730951, -9.667293313452913e-17}};

// Could move this into plink2_highprec.
dd_real ddr_log1pmx(dd_real x_ddr) {
  static const double minLog1Value = -0.79149064;
  // printf("x: %.17g\n", x_ddr.x[0]);
  if (ddr_gtd(x_ddr, 1) || (x_ddr.x[0] < minLog1Value)) {
    return ddr_sub(ddr_log1p(x_ddr), x_ddr);
  }
  const dd_real r_ddr = ddr_accurate_div(x_ddr, ddr_addd(x_ddr, 2));
  const dd_real y_ddr = ddr_sqr(r_ddr);
  if (fabs(x_ddr.x[0]) < 1e-2) {
    // |y| < (1/199)^2 < 2^{-15}; first omitted term has magnitude less than
    // ~2^{-70}|x| so this should be good enough for non-approx BASYM.
    const dd_real _ddr_2_5ths = {{4.0000000000000002220e-01, -2.2204460492503132041e-17}};
    const dd_real _ddr_2_7ths = {{2.8571428571428569843e-01, 1.5860328923216521126e-17}};
    const dd_real _ddr_2_9ths = {{2.2222222222222220989e-01, 1.2335811384723960875e-17}};
    dd_real sum_ddr = ddr_mul(_ddr_2_9ths, y_ddr);
    sum_ddr = ddr_mul(ddr_add(sum_ddr, _ddr_2_7ths), y_ddr);
    sum_ddr = ddr_mul(ddr_add(sum_ddr, _ddr_2_5ths), y_ddr);
    sum_ddr = ddr_mul(ddr_add(sum_ddr, _ddr_2_3rds), y_ddr);
    return ddr_mul(ddr_sub(sum_ddr, x_ddr), r_ddr);
  }
  // r := x/(x+2) is in [0.0944118, 1/3], y is in [0.0089135, 1/9], y/x is in
  // [0.04275, 1/9]
  // logcf term is ~1/3, so ratio between y^2 * logcfterm and x is around
  // y(y/x)(1/3) which is never much more than 1/243.  So tol_logcf=2^{-60}
  // should ensure final relative error < 2^{-67}, meeting the non-approx BASYM
  // target.
  const double tol_logcf = k2m60;
  return ddr_mul(r_ddr,
                 ddr_sub(ddr_mul(ddr_mul_pwr2(y_ddr, 2),
                                 ddr_logcf_d2(y_ddr, 3, tol_logcf)),
                         x_ddr));
}

// This term is relatively insignificant, could calculate to lower precision
// (especially when approx=True).
static inline dd_real bcorr_ddr(dd_real abmin_ddr, dd_real abmax_ddr) {
  return ddr_add(ddr_sub(ddr_stirlerr(abmax_ddr), ddr_stirlerr(ddr_add(abmin_ddr, abmax_ddr))),
                 ddr_stirlerr(abmin_ddr));
}

// must be even
CONSTI32(kBasymApproxIter, 20);

static const dd_real _ddr_half_e0_recip = {{0.443113462726379, -1.9166466249564497e-17}};  // sqrt(pi)/4

dd_real basym_approx(double a, double b, dd_real lambda_ddr) {
  const double e1 = 0.3535533905932738;  // 2^{-3/2}
  const double ln_e0 = 0.12078223763524522;

  double a0[kBasymApproxIter + 1];
  double b0[kBasymApproxIter + 1];
  double c[kBasymApproxIter + 1];
  double d[kBasymApproxIter + 1];

  // This is the dominant term if we're relatively far from the mode; it's
  // worth calculating to dd_real precision.
  // todo: try to calculate (t_ddr minus logpmf) accurately to support
  // tail-likelihood calculation
  const dd_real t_ddr = ddr_add(ddr_muld(ddr_log1pmx(ddr_divd(lambda_ddr, -a)), a),
                                ddr_muld(ddr_log1pmx(ddr_divd(lambda_ddr, b)), b));

  const dd_real f_ddr = ddr_negate(t_ddr);

  const dd_real z0_ddr = ddr_sqrt(f_ddr);
  double abmin;
  double abmax;
  if (a < b) {
    abmin = a;
    abmax = b;
  } else {
    abmin = b;
    abmax = a;
  }
  const double h = abmin / abmax;
  const double r1 = (b - a) / abmax;
  const double w0sqr = 1.0 / (abmin * (h + 1));
  const double zw0sqr = w0sqr * 2 * f_ddr.x[0];
  const double w0 = sqrt(w0sqr);
  const double r0_x2 = 2 / (h + 1);

  a0[0] = r1 * (2.0 / 3);
  d[0] = a0[0] * 0.5;
  c[0] = -d[0];
  // This is the other potential leading term.
  const dd_real initial_j0_ddr = ddr_mul(_ddr_half_e0_recip, erfcx_ddr(z0_ddr));
  double j1w = e1 * w0;
  dd_real sum_ddr = ddr_addd(initial_j0_ddr, d[0] * j1w);

  double j0w = initial_j0_ddr.x[0];

  double s = 1.0;
  const double h2 = h * h;
  double hn = 1.0;
  double e1_znm1_w0n = e1 * z0_ddr.x[0] * w0sqr * kSqrt2;
  double e1_zn_w0np1 = j1w * zw0sqr;
  for (int32_t n = 2; n <= kBasymApproxIter; n += 2) {
    hn *= h2;
    a0[n - 1] = r0_x2 * (h * hn + 1) / (n + 2);
    const int32_t np1 = n + 1;
    s += hn;
    a0[n] = r1 * 2 * s / (n + 3);

    for (int32_t i = n; i <= np1; ++i) {
      const double r = (-i - 1) * 0.5;
      b0[0] = r * a0[0];
      for (int32_t m = 2; m <= i; ++m) {
        double bsum = 0;
        for (int32_t j = 1; j < m; ++j) {
          const int32_t mmj = m - j;
          bsum = prefer_fma(prefer_fma(j, r, -mmj) * a0[j - 1], b0[mmj - 1], bsum);
        }
        b0[m - 1] = prefer_fma(r, a0[m - 1], bsum / m);
      }
      c[i - 1] = b0[i - 1] / (i + 1);

      double dsum = 0;
      for (int32_t j = 1; j < i; ++j) {
        dsum = prefer_fma(d[i - j - 1], c[j - 1], dsum);
      }
      d[i - 1] = -(dsum + c[i - 1]);
    }

    // Under the conditions we're calling this function under (where abmin is
    // close to abmax), z is proportional to sqrt(abmin) and w0 is proportional
    // to 1/sqrt(abmin).
    //
    // Original algorithm kept track of w0^n and j (a degree-n polynomial in z)
    // and multiplied by w0^n and then j; this blew up into 0 * inf = nan for
    // very large abmin.
    //
    // To fix this, we replace {w, znm1, zn, j0, j1} in the recurrence with
    // {e1_znm1_w0n, e1_zn_w0np1, j0w, j1w}.
    // Previously we had
    //   j0 = e1 * znm1 + (n - 1) * j0;
    //   j1 = e1 * zn + n * j1;
    //   znm1 = z2 * znm1;
    //   zn = z2 * zn;
    //   w *= w0;
    //   const double t0 = d[n - 1] * w * j0;
    //   w *= w0;
    //   const double t1 = d[n] * w * j1;
    j0w = prefer_fma(j0w, (n - 1) * w0sqr, e1_znm1_w0n);
    j1w = prefer_fma(j1w, n * w0sqr, e1_zn_w0np1);
    e1_znm1_w0n = zw0sqr * e1_znm1_w0n;
    e1_zn_w0np1 = zw0sqr * e1_zn_w0np1;
    const double t0 = d[n - 1] * j0w;
    const double t1 = d[n] * j1w;
    sum_ddr = ddr_add(sum_ddr, ddr_add2d(t0, t1));
    // could use wider eps when |ln_e0| is large
    // (narrowing this doesn't noticeably improve accuracy, we're limited by
    // accumulated floating-point errors in intermediate results)
    const double eps = 100 * k2m53;
    if (fabs(t0) + fabs(t1) <= eps * sum_ddr.x[0]) {
      break;
    }
  }

  return ddr_add(ddr_sub(ddr_addd(t_ddr, ln_e0), bcorr_ddr(ddr_maked(abmin), ddr_maked(abmax))), ddr_log(sum_ddr));
}

CONSTI32(kBasymIter, 26);

dd_real basym(dd_real a_ddr, dd_real b_ddr, dd_real lambda_ddr) {
  const dd_real e1_ddr = {{0.3535533905932738, -2.4168233283632284e-17}};  // 2^{-3/2}
  const dd_real ln_e0_ddr = {{0.12078223763524522, 4.1797047492946264e-18}};  // log(2/sqrt(pi))

  dd_real a0_ddr[kBasymIter + 1];
  dd_real b0_ddr[kBasymIter + 1];
  dd_real c_ddr[kBasymIter + 1];
  dd_real d_ddr[kBasymIter + 1];

  const dd_real t_ddr = ddr_add(ddr_mul(ddr_log1pmx(ddr_negate(ddr_accurate_div(lambda_ddr, a_ddr))), a_ddr),
                                ddr_mul(ddr_log1pmx(ddr_accurate_div(lambda_ddr, b_ddr)), b_ddr));

  const dd_real f_ddr = ddr_negate(t_ddr);

  const dd_real z0_ddr = ddr_sqrt(f_ddr);
  dd_real abmin_ddr;
  dd_real abmax_ddr;
  if (ddr_lt(a_ddr, b_ddr)) {
    abmin_ddr = a_ddr;
    abmax_ddr = b_ddr;
  } else {
    abmin_ddr = b_ddr;
    abmax_ddr = a_ddr;
  }
  const dd_real h_ddr = ddr_accurate_div(abmin_ddr, abmax_ddr);
  const dd_real r1_ddr = ddr_accurate_div(ddr_sub(b_ddr, a_ddr), abmax_ddr);
  const dd_real hp1_ddr = ddr_addd(h_ddr, 1);
  const dd_real w0sqr_ddr = ddr_accurate_div(ddr_maked(1), ddr_mul(hp1_ddr, abmin_ddr));
  const dd_real zw0sqr_ddr = ddr_mul(w0sqr_ddr, ddr_mul_pwr2(f_ddr, 2));
  const dd_real w0_ddr = ddr_sqrt(w0sqr_ddr);
  const dd_real r0_x2_ddr = ddr_accurate_div(ddr_maked(2), hp1_ddr);
  const dd_real r1_x2_ddr = ddr_mul_pwr2(r1_ddr, 2);

  a0_ddr[0] = ddr_mul(r1_ddr, _ddr_2_3rds);
  d_ddr[0] = ddr_mul_pwr2(a0_ddr[0], 0.5);
  c_ddr[0] = ddr_negate(d_ddr[0]);
  // this is the other potential leading term
  dd_real j0w_ddr = ddr_mul(_ddr_half_e0_recip, erfcx_ddr(z0_ddr));
  dd_real j1w_ddr = ddr_mul(e1_ddr, w0_ddr);
  dd_real sum_ddr = ddr_add(ddr_mul(d_ddr[0], j1w_ddr), j0w_ddr);

  // to explore: do we still have enough precision if we just hardcode n=2
  // iteration to use dd_reals, and then fall back to float64 afterwards?

  // const dd_real z_ddr = ddr_mul(z0_ddr, _ddr_sqrt2);
  // const dd_real z2_ddr = ddr_mul_pwr2(f_ddr, 2);
  dd_real s_ddr = ddr_maked(1.0);
  const dd_real h2_ddr = ddr_sqr(h_ddr);
  dd_real hn_ddr = ddr_maked(1.0);
  dd_real e1_znm1_w0n_ddr = ddr_mul(ddr_mul(z0_ddr, e1_ddr), ddr_mul(w0sqr_ddr, _ddr_sqrt2));
  dd_real e1_zn_w0np1_ddr = ddr_mul(j1w_ddr, zw0sqr_ddr);
  for (int32_t n = 2; n <= kBasymIter; n += 2) {
    hn_ddr = ddr_mul(hn_ddr, h2_ddr);
    a0_ddr[n - 1] = ddr_divd(ddr_mul(r0_x2_ddr, ddr_addd(ddr_mul(h_ddr, hn_ddr), 1)), n+2);
    const int32_t np1 = n+1;
    s_ddr = ddr_add(s_ddr, hn_ddr);
    a0_ddr[n] = ddr_divd(ddr_mul(r1_x2_ddr, s_ddr), n+3);

    for (int32_t i = n; i <= np1; ++i) {
      const double r = (-i - 1) * 0.5;
      b0_ddr[0] = ddr_muld(a0_ddr[0], r);
      for (int32_t m = 2; m <= i; ++m) {
        dd_real bsum_ddr = ddr_maked(0);
        for (int32_t j = 1; j < m; ++j) {
          const int32_t mmj = m - j;
          bsum_ddr = ddr_add(bsum_ddr, ddr_mul(ddr_muld(a0_ddr[j - 1], j * r - mmj), b0_ddr[mmj - 1]));
        }
        b0_ddr[m - 1] = ddr_add(ddr_muld(a0_ddr[m - 1], r), ddr_divd(bsum_ddr, m));
      }
      c_ddr[i - 1] = ddr_divd(b0_ddr[i - 1], i + 1);

      dd_real dsum_ddr = ddr_maked(0);
      for (int32_t j = 1; j < i; ++j) {
        dsum_ddr = ddr_add(dsum_ddr, ddr_mul(d_ddr[i - j - 1], c_ddr[j - 1]));
      }
      d_ddr[i - 1] = ddr_negate(ddr_add(dsum_ddr, c_ddr[i - 1]));
    }

    j0w_ddr = ddr_add(ddr_mul(j0w_ddr, ddr_muld(w0sqr_ddr, n - 1)), e1_znm1_w0n_ddr);
    j1w_ddr = ddr_add(ddr_mul(j1w_ddr, ddr_muld(w0sqr_ddr, n)), e1_zn_w0np1_ddr);
    e1_znm1_w0n_ddr = ddr_mul(zw0sqr_ddr, e1_znm1_w0n_ddr);
    e1_zn_w0np1_ddr = ddr_mul(zw0sqr_ddr, e1_zn_w0np1_ddr);
    const dd_real t0_ddr = ddr_mul(d_ddr[n - 1], j0w_ddr);
    const dd_real t1_ddr = ddr_mul(d_ddr[n], j1w_ddr);
    sum_ddr = ddr_add(sum_ddr, ddr_add(t0_ddr, t1_ddr));

    // erfcx limited to ~18-digit precision
    const double eps = k2m60;
    if (fabs(t0_ddr.x[0]) + fabs(t1_ddr.x[0]) <= eps * sum_ddr.x[0]) {
      break;
    }
  }

  // printf("%.17g %.17g %.17g %.17g\n", ln_e0_ddr.x[0], t_ddr.x[0], bcorr_ddr(a, b).x[0], sum_ddr.x[0]);
  return ddr_add(ddr_sub(ddr_add(t_ddr, ln_e0_ddr), bcorr_ddr(abmin_ddr, abmax_ddr)), ddr_log(sum_ddr));
}

// Adaptations of DiDonato and Morris's BFRAC and BASYM.  BFRAC is based on a
// continued fraction introduced in
//   Aroian LA (1941) Continued fractions for the incomplete beta function.
//   Annals of Mathematical Statistics, 12.
// BASYM is introduced in
//   DiDonato AR and Morris AH Jr (1992) Algorithm 708: Significant Digit
//   Computation of the Incomplete Beta Function Ratios.  ACM Transactions on
//   Mathematical Software, 18.
//
// For most larger cases, these expansions converge more quickly than binomial
// partial sums.  ibeta_largeab_approx() makes limited use of dd_real precision
// to address the worst precision bottlenecks, while ibeta_largeab() trades off
// speed for provably great precision.
//
// (I still have work to do in understanding the derivation and properties of
// these two expansions well enough to take a real shot at improving e.g.
// the rather similar hypergeometric cdf calculation.)
double ibeta_largeab_approx(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr, uint32_t inv, uint32_t midp_complement, uint32_t return_log) {
  // normalized always true, min(aa,bb) >= 40, max(aa,bb) >= 256
  // (usually cheaper to sum tail binomial terms directly with smaller
  // min(aa,bb); and Lanczos sum becomes less accurate)
  // now supports n up to 2^900
  // caller responsible for guaranteeing aq - bp >= 0
  //
  // * In PbinomApprox(), (a,b) is initialized to (k+1,n-k) and inv is
  //   initialized to !complement; then, if (a,p) <-> (b,q) is flipped to make
  //   aq - bp >= 0, inv is flipped.
  // * When inv is true, we return maybelog(1 - cdf(a-1)) instead of
  //   maybelog(cdf(a-1)).
  // * midp_complement is midp * (1 + complement).  When midp is true:
  //   * if inv == !complement, 0.5*pmf(a-1) is subtracted from cdf(a-1) before
  //     possible inversion.
  //   * Otherwise, 0.5*pmf(a) is added to cdf(a-1) before possible inversion.
  //
  // scipy.stats.binom.logcdf uses very similar C++ code, but is afflicted with
  // high overhead:
  //   >>> import exact_tests, scipy, timeit
  //   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
  //   0.016394874997786246
  //   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
  //   0.016445749999547843
  //   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
  //   0.016285416997561697
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

  dd_real result_ln_ddr;
  // have confirmed this is still a reasonable threshold
  const double abmin = MINV(aa, bb);
  if ((abmin < 100) || (aq_minus_bp_ddr.x[0] > abmin * 0.03)) {
    // BFRAC
    result_ln_ddr = ibeta_power_terms_d_ln(aa, bb, p_ddr, q_ddr, aq_minus_bp_ddr);
    const double result_incr = log(ibeta_continued_fraction_recip_d(aa, bb, p_ddr.x[0], q_ddr.x[0], aq_minus_bp_ddr, inv, midp_complement));
    result_ln_ddr = ddr_addd(result_ln_ddr, result_incr);
  } else {
    // BASYM
    result_ln_ddr = basym_approx(aa, bb, aq_minus_bp_ddr);
    if (midp_complement) {
      const double n = aa + bb - 1;
      const uint32_t ab_flipped = (midp_complement == inv + 1);
      const double k = aa - (!ab_flipped);
      dd_real half_pmf_ddr = ddr_add(binom_ln_prob_loader(ddr_maked(k), ddr_maked(n), p_ddr, q_ddr), _ddr_log05);
      if (ab_flipped) {
        result_ln_ddr = ddr_logspace_sub(result_ln_ddr, half_pmf_ddr);
      } else {
        result_ln_ddr = ddr_logspace_add(result_ln_ddr, half_pmf_ddr);
      }
    }
  }
  if ((!inv) && return_log) {
    return result_ln_ddr.x[0];
  }
  const dd_real result_ddr = ddr_exp(result_ln_ddr);
  if (!inv) {
    return result_ddr.x[0];
  }
  const dd_real neg_result_ddr = ddr_negate(result_ddr);
  return return_log? ddr_log1p(neg_result_ddr).x[0] : ddr_addd(neg_result_ddr, 1).x[0];
}

dd_real ibeta_continued_fraction_ddr(dd_real a_ddr, dd_real b_ddr, double n, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr) {
  // min(a,b)>=2048, n=a+b-1 can be close to DBL_MAX, aq-bp guaranteed to be at
  // least 0.03*min(a,b).
  // Due to the latter, number of continued fraction iterations should be
  // limited, but we need to be careful about overflow/underflow.
  const double cf_eps = k2m64 * (1.0 / (1 << 3));
  const dd_real aq_minus_bp_plus1_ddr = ddr_addd(aq_minus_bp_ddr, 1);
  // c = (aq-bp+1) * (a/(a+1))
  dd_real cc_ddr = ddr_mul(aq_minus_bp_plus1_ddr,
                           ddr_accurate_div(a_ddr, ddr_addd(a_ddr, 1)));
  const dd_real two_minus_p_ddr = ddr_addd(q_ddr, 1);
  dd_real ff_ddr = cc_ddr;
  dd_real dd_ddr = ddr_maked(0);
  double mm = 1.0;
  while (1) {
    const dd_real denom_ddr = ddr_addd(a_ddr, 2*mm-1);
    // m*(b-m)*p / (a+2m-1)
    // Safest to just compute (b-m)/(a+2m-1) here and handle mp later.
    const dd_real shared_frac_ddr = ddr_accurate_div(ddr_addd(b_ddr, -mm), denom_ddr);

    // ((a+m-1)/(a+2m-1)) * ((n+m) * p * shared_frac) * mp
    const dd_real cur_a_ddr = ddr_mul(ddr_mul(ddr_accurate_div(ddr_addd(a_ddr, mm-1), denom_ddr),
                                              ddr_mul(ddr_mul(ddr_add2d(n, mm), p_ddr), shared_frac_ddr)),
                                      ddr_muld(p_ddr, mm));

    // (a+m)/(a+2m+1) * (m(2-p)+(aq-bp+1)) + m + shared_frac*m*p
    const dd_real cur_b_ddr = ddr_add(ddr_mul(ddr_accurate_div(ddr_addd(a_ddr, mm), ddr_addd(a_ddr, 2*mm+1)),
                                              ddr_add(aq_minus_bp_plus1_ddr, ddr_muld(two_minus_p_ddr, mm))),
                                      ddr_addd(ddr_mul(ddr_muld(shared_frac_ddr, mm), p_ddr), mm));

    mm += 1.0;
    dd_ddr = ddr_add(ddr_mul(cur_a_ddr, dd_ddr), cur_b_ddr);
    // cur_b is at least around aq-bp, which is pretty large.
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
    // If I correctly understand what's going on here, log(delta) has
    // alternating sign and decreasing magnitude, so this should ensure less
    // than cf_eps relative error is coming from incomplete evaluation of the
    // continued fraction.  (Recall that we actually need to limit the
    // logarithm's absolute error to < ~2^{-63} to achieve the desired level of
    // accuracy when DBL_MIN < p < e^{-512} and return_log=False.)
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

double ibeta_largeab(dd_real a_ddr, dd_real b_ddr, double n, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr, uint32_t inv, uint32_t return_log) {
  dd_real result_ln_ddr;
  // caller currently ensures abmin >= 2048
  if (aq_minus_bp_ddr.x[0] > MINV(a_ddr.x[0], b_ddr.x[0]) * 0.03) {
    // (a is a synonym for k+1, b is a synonym for n-k, x is a synonym for p, y
    // is a synonym for q)
    // We want
    //   log((x^a)(y^b) / Beta(a,b))
    // = a log x + b log y + log((a+b-1)!) - log((a-1)!) - log((b-1)!)
    //
    // Since
    //   binom_ln_prob_loader(a-1, n, x, y)
    // = (a-1) log x + b log y + log((a+b-1)!) - log((a-1)!) - log(b!)
    // = (a-1) log x + b log y + log((a+b-1)!) - log((a-1)!) - log((b-1)!) - log b
    // we just need to add log(bx) to it.
    result_ln_ddr = ddr_add(binom_ln_prob_loader(ddr_addd(a_ddr, -1), ddr_maked(n), p_ddr, q_ddr),
                            ddr_log(ddr_mul(b_ddr, p_ddr)));
    // Could tighten this bound.
    if ((result_ln_ddr.x[0] < -1418.0) && (inv || (!return_log))) {
      return (return_log || (!inv))? 0.0 : 1.0;
    }
    const dd_real ff_ddr = ibeta_continued_fraction_ddr(a_ddr, b_ddr, n, p_ddr, q_ddr, aq_minus_bp_ddr);
    result_ln_ddr = ddr_sub(result_ln_ddr, ddr_log(ff_ddr));
  } else {
    result_ln_ddr = basym(a_ddr, b_ddr, aq_minus_bp_ddr);
  }
  if (!inv) {
    return return_log? result_ln_ddr.x[0] : ddr_exp(result_ln_ddr).x[0];
  }
  const dd_real neg_result_ddr = ddr_negate(ddr_exp(result_ln_ddr));
  return return_log? ddr_log1p(neg_result_ddr).x[0] : ddr_addd(neg_result_ddr, 1.0).x[0];
}

// This is a port of R 4.6.0 src/nmath/qnorm.c , which implements
//   Wichura MJ (1988) Algorithm AS 241: The Percentage Points of the Normal
//   Distribution.  Applied Statistics, 37, 477-484.
// and
//   Maechler M (2022) Asymptotic Tail Formulas For Gaussian Quantiles.
//   https://cran.r-project.org/web/packages/DPQ/vignettes/qnorm-asymp.pdf .
//
// As of this writing, the R function appears to have less-efficient polynomial
// evaluation, so it isn't yet time to make this act as a simple alias for the
// built-in R function when compiling the exactr package.
//
// Note that the R function this is derived from is GPL-2 licensed; thus this
// function cannot be used in scipy.  (Boost has a similar function which lacks
// the extended logp range.)
//
// The function name has a trailing 'D' (for double-precision) since
// plink2_stats already has a faster single-parameter QuantileToZscore()
// function which just aims for float32-like precision.
double QuantileToZscoreD(double p_or_lnp, uint32_t p_is_log) {
  int32_t sign = -1;
  if (!p_is_log) {
    if (p_or_lnp > 0.5) {
      sign = 1;
      p_or_lnp = 1.0 - p_or_lnp;
    }
  } else {
    if (p_or_lnp > -kLn2) {
      sign = 1;
      p_or_lnp = -expm1(p_or_lnp);
      p_is_log = 0;
    }
  }
  if (((!p_is_log) && (p_or_lnp >= 0.075)) || (p_is_log && (p_or_lnp >= -2.5902671654458267))) {
    const double p = p_is_log? exp(p_or_lnp) : p_or_lnp;
    const double q = 0.5 - p;
    const double r = .180625 - q * q;
    // todo: benchmark polynomial evaluation strategies
    const double numer = POLY7(r,
                               3.387132872796366608,
                               133.14166789178437745,
                               1971.5909503065514427,
                               13731.693765509461125,
                               45921.953931549871457,
                               67265.770927008700853,
                               33430.575583588128105,
                               2509.0809287301226727);
    const double denom = POLY7(r,
                               1,
                               42.313330701600911252,
                               687.1870074920579083,
                               5394.1960214247511077,
                               21213.794301586595867,
                               39307.89580009271061,
                               28729.085735721942674,
                               5226.495278852854561);
    return sign * q * numer / denom;
  }
  const double lnp = p_is_log? p_or_lnp : log(p_or_lnp);
  double r = sqrt(-lnp);
  if (r <= 5) {
    r -= 1.6;
    const double numer = POLY7(r,
                               1.42343711074968357734,
                               4.6303378461565452959,
                               5.7694972214606914055,
                               3.64784832476320460504,
                               1.27045825245236838258,
                               .24178072517745061177,
                               .0227238449892691845833,
                               7.7454501427834140764e-4);
    const double denom = POLY7(r,
                               1,
                               2.05319162663775882187,
                               1.6763848301838038494,
                               .68976733498510000455,
                               .14810397642748007459,
                               .0151986665636164571966,
                               5.475938084995344946e-4,
                               1.05075007164441684324e-9);
    return sign * numer / denom;
  }
  if (r <= 27) {
    r -= 5;
    const double numer = POLY7(r,
                               6.6579046435011037772,
                               5.4637849111641143699,
                               1.7848265399172913358,
                               .29656057182850489123,
                               .026532189526576123093,
                               .0012426609473880784386,
                               2.71155556874348757815e-5,
                               2.01033439929228813265e-7);
    const double denom = POLY7(r,
                               1,
                               .59983220655588793769,
                               .13692988092273580531,
                               .0148753612908506148525,
                               7.868691311456132591e-4,
                               1.8463183175100546818e-5,
                               1.4215117583164458887e-7,
                               2.04426310338993978564e-15);
    return sign * numer / denom;
  }
  if (r >= 6.4e8) {
    return sign * r * kSqrt2;
  }
  const double s2 = -ldexp(lnp, 1);
  double x2 = s2 - log(k2Pi * s2);
  if (r < 36000) {
    x2 = s2 - log(k2Pi * x2) - 2.0/(2.0+x2);
    if (r < 840) {
      x2 = s2 - log(k2Pi * x2) + 2*log1p(-(1 - 1/(4+x2))/(2+x2));
      if (r < 109) {
        x2 = s2 - log(k2Pi * x2) + 2*log1p(-(1 - (1 - 5/(6+x2))/(4+x2))/(2+x2));
        if (r < 55) {
          x2 = s2 - log(k2Pi * x2) + 2*log1p(-(1 - (1 - (5 - 9/(8+x2))/(6+x2))/(4+x2))/(2+x2));
        }
      }
    }
  }
  return sign * sqrt(x2);
}

#ifdef __cplusplus
}
#endif
