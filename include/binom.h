#ifndef __BINOM_H__
#define __BINOM_H__

// Binomial Exact Test library, copyright (C) 2013-2026 Christopher Chang.
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

double LnBinomCoeff(int64_t n, int64_t k);

double LnBinomMass(int64_t k, int64_t n, double p);

BoolErr BinomLnP(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, int32_t midp, double* resultp);

double BinomOneSidedLnP(int64_t obs_succ, int64_t obs_tot, double succ_odds_ratio, uint32_t succ_is_greater_alt, int32_t midp);

#ifdef __cplusplus
}
#endif

#endif  // __BINOM_H__
