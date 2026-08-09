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

#include "hypergeom_detail.h"
#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

dd_real hypergeom_ln_prob_internal(int64_t m11, int64_t m12, int64_t m21, int64_t m22) {
  const int64_t m1x = m11 + m12;
  const int64_t m2x = m21 + m22;
  const int64_t mxx = m1x + m2x;
  if (!use_tdr_for_hypergeom_lnprob(mxx)) {
    dd_real ddrs[8];
    ddrs[0] = ddr_lfact(m1x);
    ddrs[1] = ddr_lfact(m2x);
    ddrs[2] = ddr_lfact(m11 + m21);
    ddrs[3] = ddr_lfact(m12 + m22);
    ddrs[4] = ddr_negate(ddr_lfact(m11));
    ddrs[5] = ddr_negate(ddr_lfact(m12));
    ddrs[6] = ddr_negate(ddr_lfact(m21));
    ddrs[7] = ddr_negate(ddr_lfact(m22));
    return ddr_sub(ddr_sort_and_add(8, ddrs), ddr_lfact(mxx));
  } else {
    td_real tdrs[8];
    tdrs[0] = tdr_lfact(m1x);
    tdrs[1] = tdr_lfact(m2x);
    tdrs[2] = tdr_lfact(m11 + m21);
    tdrs[3] = tdr_lfact(m12 + m22);
    tdrs[4] = tdr_negate(tdr_lfact(m11));
    tdrs[5] = tdr_negate(tdr_lfact(m12));
    tdrs[6] = tdr_negate(tdr_lfact(m21));
    tdrs[7] = tdr_negate(tdr_lfact(m22));
    return ddr_make_td(tdr_sub(tdr_sort_and_add(8, tdrs), tdr_lfact(mxx)));
  }
}

void HypergeomMassMultiKPrecomp(int64_t mxx, int64_t m1x, int64_t mx1, td_real* lfact_m1x_tdr_ptr, td_real* lfact_m2x_tdr_ptr, td_real* lfact_mx1_tdr_ptr, td_real* lfact_mx2_tdr_ptr, td_real* lfact_mxx_tdr_ptr) {
  const int64_t m2x = mxx - m1x;
  const int64_t mx2 = mxx - mx1;
  if (!use_tdr_for_hypergeom_lnprob(mxx)) {
    *lfact_m1x_tdr_ptr = tdr_make_dd(ddr_lfact(m1x));
    *lfact_m2x_tdr_ptr = tdr_make_dd(ddr_lfact(m2x));
    *lfact_mx1_tdr_ptr = tdr_make_dd(ddr_lfact(mx1));
    *lfact_mx2_tdr_ptr = tdr_make_dd(ddr_lfact(mx2));
    *lfact_mxx_tdr_ptr = tdr_make_dd(ddr_lfact(mxx));
    return;
  }
  *lfact_m1x_tdr_ptr = tdr_lfact(m1x);
  *lfact_m2x_tdr_ptr = tdr_lfact(m2x);
  *lfact_mx1_tdr_ptr = tdr_lfact(mx1);
  *lfact_mx2_tdr_ptr = tdr_lfact(mx2);
  *lfact_mxx_tdr_ptr = tdr_lfact(mxx);
}

double HypergeomMassJustK(int64_t m11, int64_t mxx, int64_t m1x, int64_t mx1, const td_real lfact_m1x_tdr, const td_real lfact_m2x_tdr, const td_real lfact_mx1_tdr, const td_real lfact_mx2_tdr, const td_real lfact_mxx_tdr, uint32_t logp) {
  const int64_t mx2 = mxx - mx1;
  const int64_t m12 = m1x - m11;
  const int64_t m21 = mx1 - m11;
  const int64_t m22 = mx2 - m12;
  if (!use_tdr_for_hypergeom_lnprob(mxx)) {
    dd_real ddrs[8];
    ddrs[0] = ddr_make_td(lfact_m1x_tdr);
    ddrs[1] = ddr_make_td(lfact_m2x_tdr);
    ddrs[2] = ddr_make_td(lfact_mx1_tdr);
    ddrs[3] = ddr_make_td(lfact_mx2_tdr);
    ddrs[4] = ddr_negate(ddr_lfact(m11));
    ddrs[5] = ddr_negate(ddr_lfact(m12));
    ddrs[6] = ddr_negate(ddr_lfact(m21));
    ddrs[7] = ddr_negate(ddr_lfact(m22));
    const dd_real lnresult_ddr = ddr_sub(ddr_sort_and_add(8, ddrs), ddr_make_td(lfact_mxx_tdr));
    return logp? lnresult_ddr.x[0] : ddr_exp(lnresult_ddr).x[0];
  }
  td_real tdrs[8];
  tdrs[0] = lfact_m1x_tdr;
  tdrs[1] = lfact_m2x_tdr;
  tdrs[2] = lfact_mx1_tdr;
  tdrs[3] = lfact_mx2_tdr;
  tdrs[4] = tdr_negate(tdr_lfact(m11));
  tdrs[5] = tdr_negate(tdr_lfact(m12));
  tdrs[6] = tdr_negate(tdr_lfact(m21));
  tdrs[7] = tdr_negate(tdr_lfact(m22));
  const td_real lnresult_tdr = tdr_sub(tdr_sort_and_add(8, tdrs), lfact_mxx_tdr);
  return logp? lnresult_tdr.x[0] : tdr_exp(lnresult_tdr).x[0];
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
