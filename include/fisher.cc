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

#include "fisher.h"

#include <assert.h>
#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// *cmp_resultp is set to positive value if m11 = obs_m11 + m11_incr has higher
// likelihood than m11 = obs_m11, 0 if identical likelihood, and negative value
// if lower likelihood.
// Error is returned iff malloc fails.
// If neg_numer_ddr has not been computed yet, set its
// x[0] to DBL_MAX; it will be filled in if necessary.
BoolErr FisherCompare(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m21, uint32_t obs_m22, int32_t m11_incr, dd_real* neg_numer_ddr_ptr, intptr_t* cmp_resultp, double* dbl_ptr) {
  // Fisher 2x2 likelihood is
  //
  //   (m11+m12)! (m21+m22)! (m11+m21)! (m12+m22)!
  //   -------------------------------------------
  //     m11! m12! m21! m22! (m11+m12+m21+m22)!
  //
  // so likelihood ratio of interest is
  //
  //           obs_m11! obs_m12! obs_m21! obs_m22!
  //   ---------------------------------------------------
  //   (obs_m11+j)! (obs_m12-j)! (obs_m21-j)! (obs_m22+j)!
  //
  // where j=m11_incr.
  uint32_t numer_factorial_args[4];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m21;
  numer_factorial_args[3] = obs_m22;
  uint32_t denom_factorial_args[4];
  denom_factorial_args[0] = obs_m11 + m11_incr;
  denom_factorial_args[0] = obs_m12 - m11_incr;
  denom_factorial_args[0] = obs_m21 - m11_incr;
  denom_factorial_args[0] = obs_m22 + m11_incr;

  mp_limb_t* gmp_wkspace = nullptr;
  uintptr_t gmp_wkspace_limb_ct = 0;
  BoolErr reterr = CompareFactorialProducts(4, 0, 0, numer_factorial_args, denom_factorial_args, neg_numer_ddr_ptr, &gmp_wkspace, &gmp_wkspace_limb_ct, cmp_resultp, dbl_ptr);
  free_cond(gmp_wkspace);
  return reterr;
}

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^31.
BoolErr Fisher22LnP(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m21, uint32_t obs_m22, uint32_t midp, double* resultp) {
  // Normalize.
  if (obs_m11 > obs_m22) {
    swap_u32(&obs_m11, &obs_m22);
  }
  if (obs_m12 > obs_m21) {
    swap_u32(&obs_m12, &obs_m21);
  }
  if (S_CAST(uint64_t, obs_m11) * obs_m22 > S_CAST(uint64_t, obs_m12) * obs_m21) {
    swap_u32(&obs_m11, &obs_m12);
    swap_u32(&obs_m21, &obs_m22);
  }
  if (!midp) {
    // Fast path for p=1.
    if (S_CAST(uint64_t, obs_m11 + 1) * (obs_m22 + 1) == S_CAST(uint64_t, obs_m12) * obs_m21) {
      *resultp = 0;
      return 0;
    }
  }
  double m11 = u31tod(obs_m11);
  double m12 = u31tod(obs_m12);
  double m21 = u31tod(obs_m21);
  double m22 = u31tod(obs_m22);
  double lastp = 1;
  double tailp = 1;
  // Iterate outward to floating-point precision limit.
  while (1) {
    m12 += 1;
    m21 += 1;
    lastp *= (m11 * m22) / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
  }
  // In the common case, where we're close enough to the mode that float64
  // underflow/overflow isn't an issue, use the original algorithm: sum all
  // center relative-likelihoods, sum far-tail relative-likelihoods to
  // floating-point precision limit, return log(tailp / (tailp + centerp)).
  //
  // As with HweLnP(), we note that if we're within 172 steps of the mode and
  // the starting relative-likelihood is normalized to 1, the modal
  // relative-likelihood can be loosely bounded above by
  //   ((172^172) / 172!)^4 ~= 5.3e+292
  // which leaves enough headroom to accumulate the rest of the center-sum and
  // represent intermediate values without overflowing.
  int32_t tie_ct = 1;
  if ((S_CAST(int64_t, obs_m11) + 172) * (S_CAST(int64_t, obs_m22) + 172) >=
      (S_CAST(int64_t, obs_m12) - 172) * (S_CAST(int64_t, obs_m21) - 172)) {
    lastp = 1;
    m11 = u31tod(obs_m11);
    m12 = u31tod(obs_m12);
    m21 = u31tod(obs_m21);
    m22 = u31tod(obs_m22);
    dd_real starting_lnprob_other_component_ddr = {{DBL_MAX, 0.0}};
    double centerp = 0;
    while (1) {
      m11 += 1;
      m22 += 1;
      lastp *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      // Number of center terms is maximized with obs_m11 = 0, modal_m11 = 172,
      // other values large.
      // Since 1 + 1/2 + ... + 1/172 < 1/173 + ... + 1/53000, we're limited to
      // ~53000 terms.  Each lastp update involves 4 operations which can each
      // introduce up to 0.5 ULP relative error under the default rounding
      // mode.
      if (lastp < 1 + 53000 * 2 * k2m52) {
        if (lastp <= 1 - 53000 * 2 * k2m52) {
          tailp += lastp;
          break;
        }
        // Near-tie.  True value of lastp can be greater than, equal to, or
        // less than 1.
        const intptr_t m11_incr = S_CAST(intptr_t, m11) - obs_m11;
        intptr_t cmp_result;
        if (unlikely(FisherCompare(obs_m11, obs_m12, obs_m21, obs_m22, m11_incr, &starting_lnprob_other_component_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        if (cmp_result <= 0) {
          tailp += lastp;
          tie_ct += (cmp_result == 0);
          break;
        }
      }
      centerp += lastp;
    }
    // TODO
  }
  // TODO
  return 0;
}

#ifdef __cplusplus
}
#endif
