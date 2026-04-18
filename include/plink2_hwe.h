#ifndef __PLINK2_HWE_H__
#define __PLINK2_HWE_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

BoolErr HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp, double* resultp);

// these return 0 if close enough to Hardy-Weinberg equilibrium
BoolErr HweThresh(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh, uint32_t* out_of_eqp);

BoolErr HweThreshMidp(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh, uint32_t* out_of_eqp);

BoolErr HweThreshLnMain(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp, double ln_thresh, uint32_t* out_of_eqp);

HEADER_INLINE BoolErr HweThreshLn(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, uint32_t midp, double thresh, double ln_thresh, uint32_t* out_of_eqp) {
  // kLnNormalMin = -708.3964185...
  if (ln_thresh > -708.396) {
    if (!midp) {
      return HweThresh(obs_hets, obs_hom1, obs_hom2, thresh, out_of_eqp);
    } else {
      return HweThreshMidp(obs_hets, obs_hom1, obs_hom2, thresh, out_of_eqp);
    }
  }
  return HweThreshLnMain(obs_hets, obs_hom1, obs_hom2, midp, ln_thresh, out_of_eqp);
}

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_HWE_H__
