#ifndef __SPECIAL_FUNC_H__
#define __SPECIAL_FUNC_H__

// Copyright (C) 2013-2026 Christopher Chang.
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

#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// These functions are specialized to the range min(a,b)>=40, a+b<=2^52,
// min(p,q)>2^{-960}.  (The Boost source code is a good starting point for
// handling other ranges.)

// Returns log((p^a)(q^b) / Beta(a,b))
//       = log((p^a)(q^b)(a+b-1)! / ((a-1)!(b-1)!)).
dd_real ibeta_power_terms_d_ln(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr);

// Returns the reciprocal of the Aroian/DiDonato/Morris continued fraction
// (corresponding to a binomial tail-sum relative-likelihood) if
// midp_complement=0.  If midp_complement is nonzero, a midp-adjustment is
// performed (1 = no complement, 2 = complement).
double ibeta_continued_fraction_recip_d(double aa, double bb, double xx, double yy, dd_real ay_minus_bx_ddr, uint32_t inv, uint32_t midp_complement);

// Evaluates regularized incomplete beta function with ordinary (usually off by
// several ULPs) accuracy.
double ibeta_largeab_approx(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr, uint32_t inv, uint32_t midp_complement, uint32_t return_log);

// High-accuracy variant of ibeta_continued_fraction_recip_d().
dd_real ibeta_continued_fraction_ddr(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr);

// Evaluates regularized incomplete beta function with high (<1 ULP error)
// accuracy.
double ibeta_largeab(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real aq_minus_bp_ddr, uint32_t inv, uint32_t return_log);


double QuantileToZscoreD(double p_or_lnp, uint32_t p_is_log);

#ifdef __cplusplus
}
#endif

#endif  // __SPECIAL_FUNC_H__
