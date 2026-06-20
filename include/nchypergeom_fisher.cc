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

#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// To calculate the odds-ratio confidence intervals reported by R fisher.test,
// we need moderate-accuracy mean and cdf functions for the Fisher non-central
// hypergeometric distribution.
//
// (This is similar to a subset of Agner Fog's implementation, a BSD-licensed
// copy of which is under scipy/stats/biasedurn/ in the scipy codebase.  We
// extend the range to total < 2^52.)

int64_t ApproxModeFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds) {
  const double m1d = m1;
  const double m2d = m2;
  const double nd = n;
  if (odds == 1.0) {
    // Avoid division by zero.
    return S_CAST(int64_t, (m1d + 1) * (nd + 1) / (m1d + m2d + 2));
  }
  // Solve m11 * m22 = m12 * m21 * odds.
  const double aa = 1 - odds;
  const double bb = prefer_fma(m1d + nd, odds, m2d - nd);
  const double neg_c = odds * m1d * nd;
  // If 4ac has much smaller magnitude than b^2, we'll face a lot of
  // cancellation here.
  const dd_real discrim_ddr = ddr_add(ddr_mul2d(bb, bb), ddr_mul2d(4 * aa, neg_c));
  dd_real sqrt_discrim_ddr = ddr_maked(0.0);
  if (discrim_ddr.x[0] > 0.0) {
    sqrt_discrim_ddr = ddr_sqrt(discrim_ddr);
  }
  return S_CAST(int64_t, ddr_addd(sqrt_discrim_ddr, -bb).x[0] / (2 * aa) + 0.5);
}

dd_real MeanFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds) {
  // Start from ~mode, sum outward in both directions.
  // Don't need intermediate dd_reals to reach our accuracy target.  However,
  // when n is very large, we do need to return a dd_real just to *represent*
  // the final mean with sufficient accuracy.
  const int64_t mode = ApproxModeFNCHypergeo(m1, m2, n, odds);
  double lik = 1.0;
  double m11 = mode;
  // We never need the value of m12 in isolation.  We can get away with one
  // less multiply in the main loops if we just update (m12*odds) instead of
  // m12.  This speedup usually comes with a hit to accuracy, but given that
  // we're only using this function in the context of ~float32-precision
  // root-finding, that tradeoff is acceptable.
  double m12_odds = (m1 - mode) * odds;
  double m21 = n - mode;
  double m22 = m2 - m21;

  // Iterate rightward until convergence.
  double rnumer = 0.0;
  double rdenom = 1.0;
  // Numerator is sum((m11-mode)*p) instead of sum(m11*p), then we add mode
  // back at the end.
  double m11_minus_mode = 0;
  while (1) {
    m11 += 1;
    m22 += 1;
    lik *= (m12_odds * m21) / (m11 * m22);
    m12_odds -= odds;
    m21 -= 1;
    // rnumer converges more slowly than rdenom, so we only need to check the
    // former for convergence.
    m11_minus_mode += 1;
    const double preadd = rnumer;
    rnumer = prefer_fma(lik, m11_minus_mode, rnumer);
    // Since m12_odds can become slightly inaccurate, we're not guaranteed to
    // exit when m12_odds is supposed to hit zero.  Ensure we exit when
    // m12_odds is negative (or m21 hits zero).
    if (rnumer <= preadd) {
      rnumer = preadd;
      break;
    }
    rdenom += lik;
  }
  // Jump back to mode, and then iterate leftward until left-sums converge.
  lik = 1.0;
  m11 = mode;
  m12_odds = (m1 - mode) * odds;
  m21 = n - mode;
  m22 = m2 - m21;
  double lnumer = 0.0;
  double ldenom = 0.0;
  m11_minus_mode = 0.0;
  while (1) {
    m12_odds += odds;
    m21 += 1;
    lik *= (m11 * m22) / (m12_odds * m21);
    m11 -= 1;
    m22 -= 1;
    // ldenom converges more slowly than lnumer.
    const double preadd = ldenom;
    ldenom += lik;
    if (ldenom == preadd) {
      break;
    }
    m11_minus_mode -= 1;
    lnumer = prefer_fma(lik, m11_minus_mode, lnumer);
  }
  return ddr_add2d((lnumer + rnumer) / (ldenom + rdenom), mode);
}

double VarianceFNCHypergeoFromMean(int64_t m1, int64_t m2, int64_t n, double odds, dd_real mean_ddr) {
  // Harkness, WL (1965) Properties of the Extended Hypergeometric
  // Distribution.  Annals of Mathematical Statistics, 36.
  //
  //   (1-odds) * variance = m1*n*odds - (total - (m1+n)*(1-odds)) * mean - (1 - odds) * mean^2
  //   variance = (m1*n*odds - total*mean)/(1-odds) + mean * (m1 + n - mean)
  const double total = m1 + m2;
  const double m1d = m1;
  const double nd = n;
  if (odds == 1.0) {
    // Avoid division by zero.
    if (total < 2) {
      return 0.0;
    }
    return m1d * nd * (total - m1d) * (total - nd) / (total * total * (total - 1));
  }
  // Catastrophic cancellation possible here (e.g. m1=m2=n=2^26, odds huge,
  // mean ~= 2^26 - 1).
  if (total < (1LL << 39)) {
    const dd_real first_term_ddr = ddr_accurate_div(ddr_sub(ddr_muld(ddr_mul2d(m1d, nd), odds), ddr_muld(mean_ddr, total)), ddr_add2d(1, -odds));
    const dd_real second_term_ddr = ddr_mul(ddr_subd(mean_ddr, m1d + nd), mean_ddr);
    return ddr_sub(first_term_ddr, second_term_ddr).x[0];
  }
  const td_real mean_tdr = tdr_make_dd(mean_ddr);
  const td_real first_term_tdr = tdr_accurate_div(tdr_sub(tdr_muld(tdr_make_dd(ddr_mul2d(m1d, nd)), odds), tdr_muld(mean_tdr, total)), tdr_make_dd(ddr_add2d(1, -odds)));
  const td_real second_term_tdr = tdr_mul(tdr_subd(mean_tdr, m1d + nd), mean_tdr);
  return tdr_sub(first_term_tdr, second_term_tdr).x[0];
}

void P_FNCHypergeoTwoOdds(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double odds1, double odds2, double* result1p, double* result2p) {
  // Assumes odds1 <= odds2, and they're close to each other.
  //
  // We're only using this in a root-finding context.  Calculating this for two
  // nearby odds-ratios at a time is substantially less than twice as expensive
  // as calculating it for a single odds-ratio.
  // (obvious todo: implement a faster way to estimate derivative w.r.t. odds.)
  double lik1 = k2m52;  // avoid premature overflow
  double lik2 = k2m52;
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;

  // Iterate rightward until convergence.
  double right_sum1 = 0.0;
  double right_sum2 = 0.0;
  while (1) {
    m11 += 1;
    m22 += 1;
    const double shared_mult = (m12 * m21) / (m11 * m22);
    lik1 *= shared_mult * odds1;
    lik2 *= shared_mult * odds2;
    m12 -= 1;
    m21 -= 1;
    // right_sum2 converges more slowly than right_sum1...
    const double preadd = right_sum2;
    right_sum2 += lik2;
    if (right_sum2 == preadd) {
      break;
    }
    right_sum1 += lik1;
  }
  if (right_sum2 > DBL_MAX) {
    // ...unless it blows up.
    while (1) {
      const double preadd = right_sum1;
      right_sum1 += lik1;
      if (right_sum1 == preadd) {
        break;
      }
      m11 += 1;
      m22 += 1;
      lik1 *= odds1 * ((m12 * m21) / (m11 * m22));
      m12 -= 1;
      m21 -= 1;
    }
    if (right_sum1 > DBL_MAX) {
      *result1p = 0;
      *result2p = 0;
    }
  }
  // Jump back to mode, and then iterate leftward until left-sums converge.
  lik1 = k2m52;
  lik2 = k2m52;
  m11 = obs_m11;
  m12 = obs_m12;
  m21 = obs_m21;
  m22 = obs_m22;
  const double inv_odds1 = 1.0 / odds1;
  const double inv_odds2 = 1.0 / odds2;
  double left_sum1 = lik1;
  double left_sum2 = lik2;
  while (1) {
    m12 += 1;
    m21 += 1;
    const double shared_mult = (m11 * m22) / (m12 * m21);
    lik1 *= shared_mult * inv_odds1;
    lik2 *= shared_mult * inv_odds2;
    m11 -= 1;
    m22 -= 1;
    // left_sum1 converges more slowly than left_sum2.
    const double preadd = left_sum1;
    left_sum1 += lik1;
    if (left_sum1 == preadd) {
      break;
    }
    left_sum2 += lik2;
  }
  if (left_sum1 > DBL_MAX) {
    *result1p = 1;
  } else {
    *result1p = left_sum1 / (left_sum1 + right_sum1);
  }
  if (left_sum2 > DBL_MAX) {
    // function is unimodal, not possible for both left and right sums to
    // overflow.
    *result2p = 1;
  } else {
    *result2p = left_sum2 / (left_sum2 + right_sum2);
  }
}

#ifdef __cplusplus
}
#endif
