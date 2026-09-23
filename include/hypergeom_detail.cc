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
          m11_ddr = ddr_addd(m11_ddr, 1);
          m22 += 1;
          lik_ddr = ddr_divd(ddr_mul(ddr_accurate_div(ddr_mul(lik_ddr, m12_ddr), m11_ddr), m21_ddr), m22);
          x_ddr = ddr_add(x_ddr, lik_ddr);
          if (lik_ddr.x[0] <= lik_stop) {
            break;
          }
          m12_ddr = ddr_subd(m12_ddr, 1);
          m21_ddr = ddr_subd(m21_ddr, 1);
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
