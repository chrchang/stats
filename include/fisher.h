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

#ifdef __cplusplus
namespace plink2 {
#endif

BoolErr Fisher22TwoSidedP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, uint32_t logp, double* resultp);

double Fisher22OneSidedLnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp);

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
