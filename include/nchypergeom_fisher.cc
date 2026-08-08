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
#include "nchypergeom_fisher.h"

#ifdef __cplusplus
namespace plink2 {
#endif

dd_real nchypergeom_ln_prob_ratio(int64_t obs_m11, int64_t modal_m11, int64_t m1x, int64_t m2x, int64_t mx1, double odds) {
  if (!use_tdr_for_nchypergeom_lnprob(m1x + m2x)) {
    dd_real ddrs[9];
    ddrs[0] = ddr_lfact(modal_m11);
    ddrs[1] = ddr_lfact(m1x - modal_m11);
    int64_t m21 = mx1 - modal_m11;
    ddrs[2] = ddr_lfact(m21);
    ddrs[3] = ddr_lfact(m2x - m21);
    ddrs[4] = ddr_negate(ddr_lfact(obs_m11));
    ddrs[5] = ddr_negate(ddr_lfact(m1x - obs_m11));
    m21 = mx1 - obs_m11;
    ddrs[6] = ddr_negate(ddr_lfact(m21));
    ddrs[7] = ddr_negate(ddr_lfact(m2x - m21));
    ddrs[8] = ddr_muld(ddr_log(ddr_maked(odds)), obs_m11 - modal_m11);
    return ddr_sort_and_add(9, ddrs);
  }
  td_real tdrs[9];
  tdrs[0] = tdr_lfact(modal_m11);
  tdrs[1] = tdr_lfact(m1x - modal_m11);
  int64_t m21 = mx1 - modal_m11;
  tdrs[2] = tdr_lfact(m21);
  tdrs[3] = tdr_lfact(m2x - m21);
  tdrs[4] = tdr_negate(tdr_lfact(obs_m11));
  tdrs[5] = tdr_negate(tdr_lfact(m1x - obs_m11));
  m21 = mx1 - obs_m11;
  tdrs[6] = tdr_negate(tdr_lfact(m21));
  tdrs[7] = tdr_negate(tdr_lfact(m2x - m21));
  tdrs[8] = tdr_muld(tdr_log(tdr_make1(odds)), obs_m11 - modal_m11);
  return ddr_make_td(tdr_sort_and_add(9, tdrs));
}

intptr_t FNCHypergeoCompare(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, td_real odds_tdr, int64_t m22_incr, td_real* starting_lnprobv_tdr_ptr, td_real* ln_odds_ratio_tdr_ptr, double* dbl_ptr) {
  // Likelihood ratio of interest is
  //
  //        obs_m11! obs_m12! obs_m21! obs_m22! odds^j
  //   ---------------------------------------------------
  //   (obs_m11+j)! (obs_m12-j)! (obs_m21-j)! (obs_m22+j)!
  //
  // where j=m22_incr.
  //
  // starting_lnprobv_tdr is assumed to use odds^0, not odds^obs_m22.
  uint64_t numer_factorial_args[4];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m21;
  numer_factorial_args[3] = obs_m22;
  uint64_t denom_factorial_args[4];
  denom_factorial_args[0] = obs_m11 + m22_incr;
  denom_factorial_args[1] = obs_m12 - m22_incr;
  denom_factorial_args[2] = obs_m21 - m22_incr;
  denom_factorial_args[3] = obs_m22 + m22_incr;
  return CompareFactorialProducts(4, odds_tdr, m22_incr, 0, numer_factorial_args, denom_factorial_args, starting_lnprobv_tdr_ptr, ln_odds_ratio_tdr_ptr, dbl_ptr);
}

// To calculate the odds-ratio confidence intervals reported by R fisher.test,
// we need moderate-accuracy mean and cdf functions for the Fisher non-central
// hypergeometric distribution.
//
// (This is similar to a subset of Agner Fog's implementation, a BSD-licensed
// copy of which is under scipy/stats/biasedurn/ in the scipy codebase.  We
// extend the range to total < 2^52.)

// Could be slightly off in either direction due to floating-point error.
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

// Selects larger mode if tied.
int64_t ModeFNCHypergeo(int64_t m1x, int64_t m2x, int64_t mx1, double odds) {
  double m11 = ApproxModeFNCHypergeo(m1x, m2x, mx1, odds);
  double m12 = m1x - m11;
  double m21 = mx1 - m11;
  double m22 = m2x - m21;
  if ((m11 + 1) * (m22 + 1) <= m12 * odds * m21) {
    m11 += 1;
    m22 += 1;
    do {
      m12 -= 1;
      m21 -= 1;
      m11 += 1;
      m22 += 1;
    } while (m11 * m22 <= m12 * odds * m21);
    return m11 - 1;
  }
  while (1) {
    m12 += 1;
    m21 += 1;
    if (m11 * m22 <= m12 * odds * m21) {
      return m11;
    }
    m11 -= 1;
    m22 -= 1;
  }
}

double CentralInvProbFNCHypergeo(double obs_m11d, double obs_m12d, double obs_m21d, double obs_m22d, double odds) {
  double m11 = obs_m11d;
  double m12 = obs_m12d;
  double m21 = obs_m21d;
  double m22 = obs_m22d;
  double lik = 1;
  double left_sum = 1;
  while (1) {
    m12 += 1;
    m21 += 1;
    lik *= m11 * m22 / (m12 * m21 * odds);
    m11 -= 1;
    m22 -= 1;
    const double preadd = left_sum;
    left_sum += lik;
    if (left_sum == preadd) {
      break;
    }
  }
  m11 = obs_m11d;
  m12 = obs_m12d;
  m21 = obs_m21d;
  m22 = obs_m22d;
  lik = 1;
  double right_sum = 0;
  while (1) {
    m11 += 1;
    m22 += 1;
    lik *= (m12 * m21 * odds) / (m11 * m22);
    m12 -= 1;
    m21 -= 1;
    const double preadd = right_sum;
    right_sum += lik;
    if (right_sum == preadd) {
      break;
    }
  }
  return left_sum + right_sum;
}

// Separating the mean and variance calculations breaks down for very large n
// when (mean_ddr - mode) is only known to <53-bit accuracy; using td_real
// within VarianceFNCHypergeoFromMean doesn't work there.
/*
dd_real MeanFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds) {
  // Start from mode, sum outward in both directions.
  // Don't need intermediate dd_reals to reach our accuracy target.  However,
  // when n is very large, we do need to return a dd_real just to *represent*
  // the final mean with sufficient accuracy.
  const int64_t mode = ModeFNCHypergeo(m1, m2, n, odds);
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
*/

void MeanAndVarianceFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds, dd_real* mean_ddr_ptr, double* variance_ptr) {
  // Start from mode, sum outward in both directions.
  // Don't need intermediate dd_reals to reach our accuracy target.  However,
  // when n is very large, we do need to return mean as a dd_real just to
  // *represent* the final mean with sufficient accuracy.
  const int64_t mode = ModeFNCHypergeo(m1, m2, n, odds);
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
  double rnumer2 = 0.0;
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
    // rnumer2 converges more slowly than rnumer and rdenom, so we only need to
    // check the former for convergence.
    m11_minus_mode += 1;
    const double rnumer_incr = lik * m11_minus_mode;
    const double preadd = rnumer2;
    rnumer2 = prefer_fma(rnumer_incr, m11_minus_mode, rnumer2);
    // Since m12_odds can become slightly inaccurate, we're not guaranteed to
    // exit when m12_odds is supposed to hit zero.  Ensure we exit when
    // m12_odds is negative (or m21 hits zero).
    if (rnumer2 <= preadd) {
      rnumer2 = preadd;
      break;
    }
    rnumer += rnumer_incr;
    rdenom += lik;
  }
  // Jump back to mode, and then iterate leftward until left-sums converge.
  lik = 1.0;
  m11 = mode;
  m12_odds = (m1 - mode) * odds;
  m21 = n - mode;
  m22 = m2 - m21;
  double lnumer = 0.0;
  double lnumer2 = 0.0;
  double ldenom = 0.0;
  m11_minus_mode = 0.0;
  while (1) {
    m12_odds += odds;
    m21 += 1;
    lik *= (m11 * m22) / (m12_odds * m21);
    m11 -= 1;
    m22 -= 1;
    // ldenom converges more slowly than lnumer and lnumer2.
    const double preadd = ldenom;
    ldenom += lik;
    if (ldenom == preadd) {
      break;
    }
    m11_minus_mode -= 1;
    const double lnumer_incr = lik * m11_minus_mode;
    lnumer += lnumer_incr;
    lnumer2 = prefer_fma(lnumer_incr, m11_minus_mode, lnumer2);
  }
  const double numer = lnumer + rnumer;
  const double denom = ldenom + rdenom;
  const double shifted_mean = numer / denom;
  const double shifted_ssq = (lnumer2 + rnumer2) / denom;
  *mean_ddr_ptr = ddr_add2d(shifted_mean, mode);
  *variance_ptr = MAXV(0, shifted_ssq - shifted_mean * shifted_mean);
}

double P_FNCHypergeo(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double odds, uint32_t m11_is_greater_alt, int32_t midp, uint32_t logp) {
  // Very similar to PhyperApprox().  (We assume caller will use PhyperApprox()
  // directly when odds == 1.)

  // Normalize.
  if (obs_m11 < obs_m22) {
    swap_i64(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i64(&obs_m12, &obs_m21);
  }
  // Note that 'm11_is_greater_alt' isn't equivalent to 'complement': (if midp
  // is false) it includes the starting table while 'complement' excludes it.
  // Might want to change the interface.
  if (m11_is_greater_alt) {
    swap_i64(&obs_m11, &obs_m12);
    swap_i64(&obs_m21, &obs_m22);
    odds = 1.0 / odds;
  }
  const int64_t mode = ModeFNCHypergeo(obs_m11 + obs_m12, obs_m21 + obs_m22, obs_m11 + obs_m21, odds);
  double m11 = obs_m11;
  double m12_odds = obs_m12 * odds;
  double m21 = obs_m21;
  double m22 = obs_m22;
  if (obs_m11 >= mode) {
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know the p-value > 1 - logp? DBL_MIN : 2^{-54}
    // (at which point we just return log(1) or 1; in the logp case, don't want
    // to risk imposing a surprising denormal-handling performance penalty for
    // no good reason), or remaining left likelihoods are smaller than the
    // prevision limit.
    const double first_right_mult = m12_odds * m21 / ((m11 + 1) * (m22 + 1));
    // r + r^2 + ... = r / (1-r)
    const double right_upper_bound = 0.5 * midp + first_right_mult / (1 - first_right_mult);
    if (right_upper_bound == 0.0) {
      return logp? 0 : 1;
    }

    // Scale our starting likelihood so that we overflow to INFINITY when we'd
    // want to early-exit and return log(1) or 1; this saves us a comparison in
    // the loop.
    const double start_lik = (DBL_MAX * (logp? DBL_MIN : k2m54)) / right_upper_bound;
    double lik = start_lik;
    double left_sum = start_lik;
    while (1) {
      m12_odds += odds;
      m21 += 1;
      lik *= m11 * m22 / (m12_odds * m21);
      m11 -= 1;
      m22 -= 1;
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
    m11 = obs_m11 + 1.0;
    m12_odds = (obs_m12 - 1) * odds;
    m21 = obs_m21 - 1;
    m22 = obs_m22 + 1;
    lik = right_sum;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= m12_odds * m21 / (m11 * m22);
      m12_odds -= odds;
      m21 -= 1;
      const double preadd = right_sum;
      right_sum += lik;
      if (right_sum <= preadd) {
        right_sum = preadd;
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
  if (!logp) {
    // Early-exit is still possible in this case.
    // Computation is mostly symmetric to the right-tail case above.
    const double first_left_mult = m11 * m22 / ((m12_odds + odds) * (m21 + 1));
    // 1 + r + r^2 + ... = 1/(1-r)
    const double left_upper_bound = 1 / (1 - first_left_mult) - 0.5 * midp;
    const double start_lik = (DBL_MAX * DBL_MIN) / left_upper_bound;
    double lik = start_lik;
    double right_sum = 0;
    while (1) {
      m11 += 1;
      m22 += 1;
      lik *= m12_odds * m21 / (m11 * m22);
      m12_odds -= odds;
      m21 -= 1;
      const double preadd = right_sum;
      right_sum += lik;
      if (right_sum <= preadd) {
        right_sum = preadd;
        break;
      }
    }
    if (right_sum == INFINITY_D) {
      return 0;
    }

    double left_sum = start_lik;
    m11 = obs_m11;
    m12_odds = obs_m12 * odds;
    m21 = obs_m21;
    m22 = obs_m22;
    lik = start_lik;
    while (1) {
      m12_odds += odds;
      m21 += 1;
      lik *= m11 * m22 / (m12_odds * m21);
      m11 -= 1;
      m22 -= 1;
      const double preadd = left_sum;
      left_sum += lik;
      if (left_sum == preadd) {
        break;
      }
    }
    const double midp_numer = -0.5 * midp * start_lik;
    const double denom = left_sum + right_sum;
    return (left_sum + midp_numer) / denom;
  }
  // We're to the left of the mode, and are responsible for tiny p-values.
  // Can't avoid evaluating right_sum to precision limit, starting from the
  // mode.
  // There are still two subcases:
  // - We naturally bump into the start of left_sum while evaluating right_sum.
  //   This is similar to the previous case, just with no overflow/early-exit
  //   possibility.
  // - Evaluation of right_sum (as a multiple of pmf(mode)) completes without
  //   hitting the observed contingency table.  Then we compute log(pmf(obs) /
  //   pmf(mode)), and evaluate left_sum as a multiple of pmf(obs).
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx1 = obs_m11 + obs_m21;
  m11 = mode;
  m12_odds = (m1x - m11) * odds;
  m21 = mx1 - m11;
  m22 = m2x - m21;
  double lik = 1;
  double right_sum = 1;
  while (1) {
    m11 += 1;
    m22 += 1;
    lik *= m12_odds * m21 / (m11 * m22);
    m12_odds -= odds;
    m21 -= 1;
    const double preadd = right_sum;
    right_sum += lik;
    if (right_sum <= preadd) {
      right_sum = preadd;
      break;
    }
  }
  const double obs_m11_d = obs_m11;
  m11 = mode;
  m12_odds = (m1x - m11) * odds;
  m21 = mx1 - m11;
  m22 = m2x - m21;
  lik = 1;
  while (1) {
    m12_odds += odds;
    m21 += 1;
    lik *= m11 * m22 / (m12_odds * m21);
    m11 -= 1;
    m22 -= 1;
    if (m11 == obs_m11_d) {
      double left_sum = lik;
      while (1) {
        m12_odds += odds;
        m21 += 1;
        lik *= m11 * m22 / (m12_odds * m21);
        m11 -= 1;
        m22 -= 1;
        const double preadd = left_sum;
        left_sum += lik;
        if (left_sum == preadd) {
          break;
        }
      }
      const double midp_numer = -0.5 * midp;
      const double denom = left_sum + right_sum;
      return log((left_sum + midp_numer) / denom);
    }
    const double preadd = right_sum;
    right_sum += lik;
    if (right_sum == preadd) {
      break;
    }
  }
  const dd_real lnprob_ratio_ddr = nchypergeom_ln_prob_ratio(obs_m11, mode, m1x, m2x, mx1, odds);
  m11 = obs_m11;
  m12_odds = (m1x - m11) * odds;
  m21 = mx1 - m11;
  m22 = m2x - m21;
  lik = 1;
  double left_sum = 1 - 0.5 * midp;
  while (1) {
    m12_odds += odds;
    m21 += 1;
    lik *= m11 * m22 / (m12_odds * m21);
    m11 -= 1;
    m22 -= 1;
    const double preadd = left_sum;
    left_sum += lik;
    if (left_sum == preadd) {
      break;
    }
  }
  return join_log_and_nonlog(lnprob_ratio_ddr, left_sum / right_sum, 1);
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
