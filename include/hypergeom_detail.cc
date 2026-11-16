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

#include <math.h>

#include "hypergeom_detail.h"
#include "plink2_float.h"
#include "plink2_highprec.h"
#include "special_func.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Should always have <1 ULP error; and ddr_exp(result) also has <1 ULP error
// when it isn't < DBL_MIN.
dd_real hypergeom_ln_prob_loader(dd_real m11_ddr, dd_real m12_ddr, dd_real m21_ddr, dd_real m22_ddr) {
  // Catherine Loader's algorithm, see R src/nmath/dhyper.c .  We subtract a
  // binomial log-probability from the sum of two others, where the magnitude
  // of the subtracted value is bounded above by 0.5 * log(mxx) < 312.
  // binom_ln_prob_loader() is accurate to ~90+ bits, so absolute error is
  // limited to ~2^{-80}, and cancellation is not a concern unless the
  // log-probability is very close but unequal to 0.
  //
  // The very-close-but-unequal-to-0 case requires exactly one cell to be 0
  // (WLOG let that be m22), and m12*m21 << m11.  The cancellation problem in
  // this scenario can be seen in R 4.6.1: stats::dhyper(0, 1, 2**52 - 2, 1,
  // log=TRUE) is ~1.5x the true value, and stats::dhyper(0, 1, 2**55, 1,
  // log=TRUE) is 0 when the true value is around -2^{-55}.  Fortunately, this
  // case is straightforward to handle by other means.

  // Normalize: m11 >= m22, m12 >= m21, m21 >= m22.
  if (ddr_lt(m11_ddr, m22_ddr)) {
    swap_ddr(&m11_ddr, &m22_ddr);
  }
  if (ddr_lt(m12_ddr, m21_ddr)) {
    swap_ddr(&m12_ddr, &m21_ddr);
  }
  if (ddr_lt(m21_ddr, m22_ddr)) {
    swap_ddr(&m11_ddr, &m12_ddr);
    swap_ddr(&m21_ddr, &m22_ddr);
  }
  if (ddr_is_zero(m22_ddr)) {
    // avoid overflow/underflow when possible
    if (ddr_is_zero(m21_ddr)) {
      return ddr_maked(0);
    }
    const double approx_first_mult = (m12_ddr.x[0] / (m11_ddr.x[0] + 1)) * m21_ddr.x[0];
    if (approx_first_mult < (1.0 / (1 << 13))) {
      // We've hit the edge case, adjacent table has probability < 2^{-13} the
      // starting contingency table.  (We're targeting an epsilon of 2^{-67}
      // for non-approximate functions, 2^{-80} / 2^{-67} = 2^{-13}.)
      //
      // Compute x := (pmf(1) + pmf(2) + ...)/pmf(0).  Then, pmf(0) = 1/(1+x) =
      // 1 - x/(1+x), so we can evaluate log(pmf(0)) as log1p(-x/(1+x)).
      if (approx_first_mult < DBL_MIN * (1 + 3 * k2m52)) {
        return ddr_maked(0);
      }
      m11_ddr = ddr_addd(m11_ddr, 1);
      dd_real x_ddr = ddr_mul(ddr_accurate_div(m12_ddr, m11_ddr), m21_ddr);
      if (x_ddr.x[0] > (k2m64 / 16)) {
        dd_real lik_ddr = x_ddr;
        const double lik_stop = x_ddr.x[0] * (k2m64 / 16);
        double m22 = 1;
        while (1) {
          // bugfix (16 Nov 2026): m12 and m21 must be decremented before the
          // multiplication, not after; pmf(j+1)/pmf(j) =
          // (m12-j)(m21-j) / ((m11+j+1)(j+1)).
          m11_ddr = ddr_addd(m11_ddr, 1);
          m12_ddr = ddr_subd(m12_ddr, 1);
          m21_ddr = ddr_subd(m21_ddr, 1);
          m22 += 1;
          lik_ddr = ddr_divd(ddr_mul(ddr_accurate_div(ddr_mul(lik_ddr, m12_ddr), m11_ddr), m21_ddr), m22);
          x_ddr = ddr_add(x_ddr, lik_ddr);
          if (lik_ddr.x[0] <= lik_stop) {
            break;
          }
        }
      }
      return ddr_log1p(ddr_negate(ddr_accurate_div(x_ddr, ddr_addd(x_ddr, 1))));
    }
  }
  const dd_real m1x_ddr = ddr_add(m11_ddr, m12_ddr);
  const dd_real m2x_ddr = ddr_add(m21_ddr, m22_ddr);
  const dd_real mxx_ddr = ddr_add(m1x_ddr, m2x_ddr);
  const dd_real q_ddr = ddr_accurate_div(m1x_ddr, mxx_ddr);
  const dd_real p_ddr = ddr_negate(ddr_addd(q_ddr, -1));
  const dd_real p1_ddr = binom_ln_prob_loader(m21_ddr, ddr_add(m11_ddr, m21_ddr), p_ddr, q_ddr);
  const dd_real p2_ddr = binom_ln_prob_loader(m22_ddr, ddr_add(m12_ddr, m22_ddr), p_ddr, q_ddr);
  const dd_real p3_ddr = binom_ln_prob_loader(m2x_ddr, mxx_ddr, p_ddr, q_ddr);
  return ddr_sub(ddr_add(p1_ddr, p2_ddr), p3_ddr);
}

static const double kTwoPi = 6.283185307179586;

// ln(n!) - ln(sqrt(2*pi*n)*(n/e)^n), for positive integer n.
static double stirlerr_dbl(double nn) {
  if (nn < 256) {
    return ddr_stirlerr(ddr_maked(nn)).x[0];
  }
  // For n >= 256, the first omitted term (n^{-9}/1188) is < 2^{-75}, while the
  // return value is > 2^{-12}.
  const double invn = 1.0 / nn;
  const double invn2 = invn * invn;
  return ((((-1.0 / 1680) * invn2 + (1.0 / 1260)) * invn2 - (1.0 / 360)) * invn2 + (1.0 / 12)) * invn;
}

// x log(x/M) + M - x, where M = m_hi + m_lo is represented exactly as an
// unevaluated sum.  x must be positive.
//
// Keeping M at double-double precision matters: the derivative of this
// function w.r.t. M is (M-x)/M, so rounding M to double precision would add
// ~2^{-53}|M-x| absolute error, which is large when M is large and x is
// several standard deviations away.
static double loader_bd0_dbl(double xx, double m_hi, double m_lo) {
  // x - m_hi is exact when they're within a factor of 2 of each other.
  const double diff = (xx - m_hi) - m_lo;
  const double sum = xx + m_hi;
  // Loader's original algorithm switches to the direct formula when |v| >=
  // 0.1, where v := (x-M)/(x+M).  But the direct formula suffers from
  // cancellation, which costs us several bits when |v| is in [0.1, 0.5) (this
  // matters for small/unbalanced cells).  The series converges by a factor of
  // v^2 per term, so it's still cheap with |v| < 0.5.
  if (fabs(diff) >= 0.5 * sum) {
    // log(x/M) = log(x/m_hi) - log1p(m_lo/m_hi) ~= log(x/m_hi) - m_lo/m_hi
    return xx * (log(xx / m_hi) - m_lo / m_hi) - diff;
  }
  const double vv = diff / sum;
  const double v2 = vv * vv;
  double ss = diff * vv;
  double ej = 2 * xx * vv;
  for (double jj = 1; ; jj += 1) {
    ej *= v2;
    const double s1 = ss + ej / (2 * jj + 1);
    if (s1 == ss) {
      return ss;
    }
    ss = s1;
  }
}

static inline double loader_bd0_pair_dbl(double kk, double nn, double pp, double qq) {
  double np_lo;
  const double np_hi = qd_two_prod(nn, pp, &np_lo);
  double nq_lo;
  const double nq_hi = qd_two_prod(nn, qq, &nq_lo);
  return loader_bd0_dbl(kk, np_hi, np_lo) + loader_bd0_dbl(nn - kk, nq_hi, nq_lo);
}

double hypergeom_ln_prob_dbl(double m11, double m12, double m21, double m22) {
  // Same decomposition as hypergeom_ln_prob_loader(): we sum three
  // binomial log-probabilities (with a shared p, which cancels out exactly),
  // each expressed via Loader's algorithm.  To avoid cancellation, we
  // combine the three binomials' stirlerr() terms, bd0() terms, and
  // 0.5 * log(2 * pi * k * (n-k) / n) terms separately.  The latter combine
  // into a single log() call.
  //
  // Absolute error is a small multiple of 2^{-53} times
  //   1 + (sum of bd0() term magnitudes),
  // and the bd0() terms are on the order of z^2 / 2, where z is the number of
  // standard deviations between the table and the distribution mean.

  // Normalize: m11 >= m22, m12 >= m21, m21 >= m22.
  if (m11 < m22) {
    const double tmp = m11;
    m11 = m22;
    m22 = tmp;
  }
  if (m12 < m21) {
    const double tmp = m12;
    m12 = m21;
    m21 = tmp;
  }
  if (m21 < m22) {
    double tmp = m11;
    m11 = m12;
    m12 = tmp;
    tmp = m21;
    m21 = m22;
    m22 = tmp;
  }
  const double m1x = m11 + m12;
  const double m2x = m21 + m22;
  const double mx1 = m11 + m21;
  const double mx2 = m12 + m22;
  const double mxx = m1x + m2x;
  const double pp = m2x / mxx;
  const double qq = m1x / mxx;
  if (m22 == 0) {
    if (m21 == 0) {
      return 0;
    }
    const double approx_first_mult = (m12 / (m11 + 1)) * m21;
    if (approx_first_mult < (1.0 / (1 << 13))) {
      // log-probability is very close to zero; see hypergeom_ln_prob_loader().
      if (approx_first_mult < DBL_MIN * (1 + 3 * k2m52)) {
        return 0;
      }
      double xx = approx_first_mult;
      double lik = xx;
      const double lik_stop = xx * k2m54;
      double cur_m11 = m11 + 1;
      double cur_m12 = m12;
      double cur_m21 = m21;
      double cur_m22 = 1;
      while (lik > lik_stop) {
        cur_m11 += 1;
        cur_m22 += 1;
        cur_m12 -= 1;
        cur_m21 -= 1;
        lik *= (cur_m12 / cur_m11) * (cur_m21 / cur_m22);
        xx += lik;
      }
      return log1p(-xx / (1 + xx));
    }
    // Second binomial term has k=0, so it's just n log q.
    if (m11 == 0) {
      // First binomial term has k=n, so it's just n log p.
      const double stirlerr_sum = stirlerr_dbl(m1x) + stirlerr_dbl(m2x) - stirlerr_dbl(mxx);
      const double half_lf = 0.5 * log(kTwoPi * (m2x / mxx) * m1x);
      return (mx1 * log(pp) + mx2 * log(qq)) + ((stirlerr_sum + loader_bd0_pair_dbl(m2x, mxx, pp, qq)) + half_lf);
    }
    const double stirlerr_sum = (stirlerr_dbl(mx1) - stirlerr_dbl(m11) - stirlerr_dbl(m21)) - (stirlerr_dbl(mxx) - stirlerr_dbl(m1x) - stirlerr_dbl(m2x));
    const double bd0_sum = loader_bd0_pair_dbl(m21, mx1, pp, qq) - loader_bd0_pair_dbl(m2x, mxx, pp, qq);
    // (Only two binomial terms here, so the 2*pi factors cancel.)
    const double half_lf_sum = 0.5 * log(((m21 / mx1) * (m11 / m1x)) * (mxx / m2x));
    return mx2 * log(qq) + ((stirlerr_sum - bd0_sum) - half_lf_sum);
  }
  const double stirlerr_sum = (stirlerr_dbl(mx1) - stirlerr_dbl(m11) - stirlerr_dbl(m21)) + (stirlerr_dbl(mx2) - stirlerr_dbl(m12) - stirlerr_dbl(m22)) - (stirlerr_dbl(mxx) - stirlerr_dbl(m1x) - stirlerr_dbl(m2x));
  const double bd0_sum = (loader_bd0_pair_dbl(m21, mx1, pp, qq) + loader_bd0_pair_dbl(m22, mx2, pp, qq)) - loader_bd0_pair_dbl(m2x, mxx, pp, qq);
  //   0.5 * log(2 * pi * m21 * m11 / mx1)
  // + 0.5 * log(2 * pi * m22 * m12 / mx2)
  // - 0.5 * log(2 * pi * m2x * m1x / mxx)
  const double half_lf_sum = 0.5 * log(kTwoPi * ((m21 / mx1) * (m11 / m1x)) * ((m22 / mx2) * (m12 / m2x)) * mxx);
  return (stirlerr_sum - bd0_sum) - half_lf_sum;
}

intptr_t HypergeomCompare(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, int64_t m22_incr, td_real* neg_numer_tdr_ptr, double* dbl_ptr) {
  // Likelihood ratio of interest is
  //
  //           obs_m11! obs_m12! obs_m21! obs_m22!
  //   ---------------------------------------------------
  //   (obs_m11+j)! (obs_m12-j)! (obs_m21-j)! (obs_m22+j)!
  //
  // where j=m22_incr.
  //
  // Note that HWE kind of maps to this, via
  //   m11 := obs_hets*0.5
  //   m12 := obs_hom1
  //   m21 := obs_hom2
  //   m22 := (obs_hets-1)*0.5
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
  td_real ln_odds_ratio_tdr = tdr_make1(0.0);
  return CompareFactorialProducts(4, tdr_make1(1.0), 0, 0, numer_factorial_args, denom_factorial_args, neg_numer_tdr_ptr, &ln_odds_ratio_tdr, dbl_ptr);
}

#ifdef __cplusplus
}
#endif
