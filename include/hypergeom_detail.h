#ifndef __HYPERGEOM_DETAIL_H__
#define __HYPERGEOM_DETAIL_H__

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

// Support routines specific to the hypergeometric distribution.
//
// m11/m12/m21/m22 refer to cells in a 2x2 contingency table.

// Slightly smaller than use_tdr_for_binom_lnprob() threshold since the sum has
// more terms.
HEADER_INLINE uint32_t use_tdr_for_hypergeom_lnprob(int64_t obs_tot) {
  return (obs_tot >= (1LL << 35));
}

// Returns log of
//
//   (m11+m12)! (m21+m22)! (m11+m21)! (m12+m22)!
//   -------------------------------------------
//     m11! m12! m21! m22! (m11+m12+m21+m22)!
dd_real hypergeom_ln_prob_internal(int64_t m11, int64_t m12, int64_t m21, int64_t m22);

// Returns positive value if m22 = obs_m22 + m22_incr has higher probability
// than m22 = obs_m22, 0 if identical probability, and negative value if lower
// probability.
// If neg_numer_tdr has not been computed yet, set its x[0] to DBL_MAX; it will
// be filled in if necessary.
intptr_t HypergeomCompare(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, int64_t m22_incr, td_real* neg_numer_tdr_ptr, double* dbl_ptr);

HEADER_INLINE intptr_t HypergeomCompareDdr(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, int64_t m22_incr, dd_real neg_numer_ddr, double* dbl_ptr) {
  td_real neg_numer_tdr = {{neg_numer_ddr.x[0], neg_numer_ddr.x[1], DBL_MAX}};
  return HypergeomCompare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &neg_numer_tdr, dbl_ptr);
}

#ifdef __cplusplus
}
#endif

#endif  // __HYPERGEOM_DETAIL_H__
