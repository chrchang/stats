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

int64_t ApproxModeFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds);

dd_real MeanFNCHypergeo(int64_t m1, int64_t m2, int64_t n, double odds);

double VarianceFNCHypergeoFromMean(int64_t m1, int64_t m2, int64_t n, double odds, dd_real mean_ddr);

HEADER_INLINE double MeanFNCHypergeoDerivOdds(int64_t m1, int64_t m2, int64_t n, double odds, dd_real mean_ddr) {
  return VarianceFNCHypergeoFromMean(m1, m2, n, odds, mean_ddr) / odds;
}

void P_FNCHypergeoTwoOdds(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double odds1, double odds2, double* result1p, double* result2p);

#ifdef __cplusplus
}
#endif

#endif  // __NCHYPERGEOM_FISHER_H__
