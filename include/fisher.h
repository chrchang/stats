#ifndef __FISHER_H__
#define __FISHER_H__

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

#include "plink2_base.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

dd_real fisher22_ln_prob_internal(int64_t m11, int64_t m12, int64_t m21, int64_t m22);

HEADER_INLINE double HypergeomMass(int64_t m11, int64_t m12, int64_t m21, int64_t m22, uint32_t logp) {
  const dd_real ln_prob_ddr = fisher22_ln_prob_internal(m11, m12, m21, m22);
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  return ddr_exp(ln_prob_ddr).x[0];
}

BoolErr Fisher22TwoSidedP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, uint32_t logp, double* resultp);

double PhyperApprox(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp, uint32_t logp);

double Phyper(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, uint32_t logp);

int64_t Qhyper(dd_real p_or_lnp_ddr, int64_t ac, int64_t bd, int64_t ab, uint32_t logp);

HEADER_INLINE int64_t QhyperHalfUlp(dd_real p_or_lnp_ddr, int64_t ac, int64_t bd, int64_t ab, uint32_t logp) {
  if (!ddr_is_zero(p_or_lnp_ddr)) {
    const double half_ulp = 0.5 * fabs(p_or_lnp_ddr.x[0] - prev_float64(p_or_lnp_ddr.x[0]));
    p_or_lnp_ddr = ddr_subd(p_or_lnp_ddr, half_ulp);
  }
  return Qhyper(p_or_lnp_ddr, ac, bd, ab, logp);
}

BoolErr Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp, double* resultp);

// Probable todos:
// - Function implementing 2xk for k>2, recursion can just follow how
//   Fisher23LnP builds on 2x2
// - Function implementing jxk for 3<=j<=k; recursive function can take a
//   parameter specifying the column number in the bottom row to iterate over
//   (bottom-row values further to the right are locked during the function
//   call).  I expect the function would rarely be worth using over simulation
//   when j>3 (complexity is roughly O(n^{(j-1)(k-1)/2})), but O(n^2) for j=k=3
//   seems likely to be a good deal.

#ifdef __cplusplus
}
#endif

#endif  // __FISHER_H__
