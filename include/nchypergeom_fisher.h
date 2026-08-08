#ifndef __NCHYPERGEOM_FISHER_H__
#define __NCHYPERGEOM_FISHER_H__

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

// odds is assumed to be in [2^{-400}, 2^400].  (Quadratic formula in
// ApproxModeFNCHypergeo() becomes vulnerable to overflow a bit beyond that
// point.)

HEADER_INLINE uint32_t use_tdr_for_nchypergeom_lnprob(int64_t obs_tot) {
  return (obs_tot >= (1LL << 35));
}

intptr_t FNCHypergeoCompare(uint64_t obs_m11, uint64_t obs_m12, uint64_t obs_m21, uint64_t obs_m22, td_real odds_tdr, int64_t m22_incr, td_real* starting_lnprobv_tdr_ptr, td_real* ln_odds_ratio_tdr_ptr, double* dbl_ptr);

int64_t ModeFNCHypergeo(int64_t m1x, int64_t m2x, int64_t mx1, double odds);

// Just treats starting contingency table as likelihood 1, and sums leftward
// and rightward to precision limit.  Caller is responsible for not using this
// function when it might overflow.
double CentralInvProbFNCHypergeo(double obs_m11d, double obs_m12d, double obs_m21d, double obs_m22d, double odds);

// dd_real MeanFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds);

// double VarianceFNCHypergeoFromMean(int64_t m1, int64_t m2, int64_t n, double odds, dd_real mean_ddr);

void MeanAndVarianceFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds, dd_real* mean_ddr_ptr, double* variance_ptr);

double P_FNCHypergeo(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double odds, uint32_t m11_is_greater_alt, int32_t midp, uint32_t logp);

void P_FNCHypergeoTwoOdds(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double odds1, double odds2, double* result1p, double* result2p);

#ifdef __cplusplus
}
#endif

#endif  // __NCHYPERGEOM_FISHER_H__
