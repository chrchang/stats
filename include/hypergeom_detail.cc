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
