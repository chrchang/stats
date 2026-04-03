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

BoolErr Fisher22LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, uint32_t midp, double* resultp) {
  // TODO
  return 0;
}

#ifdef __cplusplus
}
#endif
