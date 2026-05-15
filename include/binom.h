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
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

double LnBinomCoeff(int64_t n, int64_t k);

double BinomMass(int64_t k, int64_t n, dd_real p_ddr, uint32_t logp);

BoolErr BinomP(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, int32_t midp, uint32_t logp, double* resultp);

double PbinomApprox(int64_t obs_k, int64_t n, dd_real p_ddr, uint32_t complement, int32_t midp, uint32_t logp);

HEADER_INLINE double BinomOneSidedP(int64_t obs_k, int64_t n, dd_real p_ddr, uint32_t succ_is_greater_alt, int32_t midp, uint32_t logp) {
  const int64_t k_decr = succ_is_greater_alt && (!midp);
  if (k_decr && (obs_k == 0)) {
    return logp? 0.0 : 1.0;
  }
  return PbinomApprox(obs_k - k_decr, n, p_ddr, succ_is_greater_alt, midp, logp);
}

double Pbinom(int64_t obs_k, int64_t n, dd_real p_ddr, uint32_t complement, uint32_t logp);

int64_t Qbinom(dd_real targetp_ddr, int64_t n, dd_real distp_ddr, uint32_t log_target);

HEADER_INLINE int64_t QbinomHalfUlp(dd_real targetp_ddr, int64_t n, dd_real distp_ddr, uint32_t log_target) {
  const double half_ulp = 0.5 * (targetp_ddr.x[0] - nextafter(targetp_ddr.x[0], 0.0));
  targetp_ddr = ddr_subd(targetp_ddr, half_ulp);
  return Qbinom(targetp_ddr, n, distp_ddr, log_target);
}

#ifdef __cplusplus
}
#endif

#endif  // __BINOM_H__
