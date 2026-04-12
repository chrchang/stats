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

#include "binom.h"

#include <assert.h>
#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Assumes exponent > 0.
void mul_by_u63_to_u31_pow(uint64_t base, uint32_t exponent, mp_limb_t* num, uint32_t* num_limb_ctp, mp_limb_t* tmppow, mp_limb_t* numcopy) {
  uint32_t num_limb_ct = *num_limb_ctp;
  uint32_t tmppow_limb_ct = 1;
  uint32_t tmppow_exp = 1;
  tmppow[0] = base;
#ifndef __LP64__
  if (base > 0xffffffffU) {
    tmppow[1] = base >> 32;
    ++tmppow_limb_ct;
  }
#endif
  while (1) {
    if (tmppow_exp & exponent) {
      memcpy(numcopy, num, num_limb_ct * sizeof(mp_limb_t));
      mp_limb_t msl;
      if (num_limb_ct >= tmppow_limb_ct) {
        msl = mpn_mul(num, numcopy, num_limb_ct, tmppow, tmppow_limb_ct);
      } else {
        msl = mpn_mul(num, tmppow, tmppow_limb_ct, numcopy, num_limb_ct);
      }
      num_limb_ct += tmppow_limb_ct - (msl == 0);
      exponent -= tmppow_exp;
      if (!exponent) {
        return;
      }
    }
    memcpy(numcopy, tmppow, tmppow_limb_ct * sizeof(mp_limb_t));
    mpn_sqr(tmppow, numcopy, tmppow_limb_ct);
    tmppow_limb_ct *= 2;
    if (tmppow[tmppow_limb_ct - 1] == 0) {
      --tmppow_limb_ct;
    }
    exponent = exponent * 2;
  }
}


// Forked from CompareFactorialProducts() in plink2_highprec; it will replace
// the original if plink2 ever needs support for a non-power-of-2 odds ratio.
//
// Preconditions:
// - numer_factorial_args[] and denom_factorial_args[] are
//   not-necessarily-sorted lists of length ffac_ct, describing a quotient of
//   factorial-products.  (If one list is longer than the other, just pad the
//   other with zeroes.)
// - odds_ratio_{numer,denom,pow} describe a number of the form
//     (odds_ratio_numer/odds_ratio_denom)^odds_ratio_pow
//   to multiply the quotient by at the end.
// - starting_lnprobv_ddr is either initialized to
//     log(odds_ratio^numer_odds_ratio_pow / (numer_factorial_args[0]! ... numer_factorial_args[ffac_ct-1]!))
//   or it has x[0] initialized to DBL_MAX to indicate that the calculation
//   hasn't happened.  In the latter case, it may be set to the former value if
//   that is needed.
// - Similarly, ln_odds_ratio_ddr is either initialized to log(odds_ratio), or
//   it has x[0] initialized to DBL_MAX.
//
// This function errors out iff memory allocation fails.
//
// Postconditions on success:
// - *cmp_resultp is set to a positive value if the fraction > 1, a negative
//   value if the fraction < 1, and zero if it's exactly 1.
// - *dbl_ptr is the double representation of the fraction, error limited to
//   1-2 ulps.
// - numer_factorial_args[] and denom_factorial_args[] are sorted in
//   nondecreasing order.
//
// This could take a precomputed ddr_lfact table as an additional pair of
// parameters, but I don't think that makes much of a difference.
BoolErr CompareFactorialProductsEx(uint32_t ffac_ct, int64_t odds_ratio_numer, int64_t odds_ratio_denom, int32_t odds_ratio_pow, int32_t numer_odds_ratio_pow, uint32_t* numer_factorial_args, uint32_t* denom_factorial_args, dd_real* starting_lnprobv_ddr_ptr, dd_real* ln_odds_ratio_ddr_ptr, mp_limb_t** gmp_wkspacep, uintptr_t* gmp_wkspace_limb_ctp, intptr_t* cmp_resultp, double* dbl_ptr) {
  // 1. Sort numer_factorial_args and denom_factorial_args.  (This has the
  //    effect of cancelling out matching terms.)
  // 2. Iterate through numer_factorial_args[] and denom_factorial_args[] just
  //    to determine bignum calculation size.
  // 3. If bignum calculation size is large enough, perform the comparison with
  //    dd_reals first, returning a result unless the log-likelihoods are
  //    within ffac_ct * 2^{-60} of each other.
  // 4. Reallocate workspace if necessary, and iterate properly through
  //    numer_factorial_args[] and denom_factorial_args[].  When
  //    numer_factorial_args[k] > denom_factorial_args[k], multiply the
  //    numerator by the corresponding falling-factorial; and vice versa.
  // 5. Multiply by power of odds_ratio.
  // 6. Perform final left-shift of numerator or denominator (still necessary
  //    since we pull out some factors of 2 while accumulating the
  //    falling-factorial products).
  // 7. Compare.
  //
  // Possible improvement: find appropriate spots to use mpn_sqr().
  STD_SORT(ffac_ct, u32cmp, numer_factorial_args);
  STD_SORT(ffac_ct, u32cmp, denom_factorial_args);

  uintptr_t numer_term_ct = 0;
  uintptr_t denom_term_ct = 0;
  uint32_t max_ffac_size = 0;
  for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
    const uint32_t numer_factorial_arg = numer_factorial_args[ffac_idx];
    const uint32_t denom_factorial_arg = denom_factorial_args[ffac_idx];
    uint32_t ffac_size;
    if (numer_factorial_arg > denom_factorial_arg) {
      ffac_size = numer_factorial_arg - denom_factorial_arg;
      numer_term_ct += ffac_size;
    } else {
      ffac_size = denom_factorial_arg - numer_factorial_arg;
      denom_term_ct += ffac_size;
    }
    if (max_ffac_size < ffac_size) {
      max_ffac_size = ffac_size;
    }
  }
  if (max_ffac_size == 0) {
    // All factorials cancel out.
    if ((odds_ratio_numer == odds_ratio_denom) || (odds_ratio_pow == 0)) {
      *dbl_ptr = 1;
      *cmp_resultp = 0;
      return 0;
    }
    *cmp_resultp = ((odds_ratio_numer > odds_ratio_denom) == (odds_ratio_pow > 0))? 1 : -1;
    // We want to guarantee *dbl_ptr isn't off by more than 3 ULPs (at least
    // when the ratio is in [0.5, 2]), so we don't use plain double arithmetic
    // past |odds_ratio_pow| == 1.
    if ((odds_ratio_pow == -1) || (odds_ratio_pow == 1)) {
      if (odds_ratio_pow == -1) {
        swap_i64(&odds_ratio_numer, &odds_ratio_denom);
      }
      // We allow odds_ratio_numer and odds_ratio_denom > 2^53, so each
      // cast-to-double can have 0.5 ULP rounding error, and then the division
      // is another potential 0.5 ULP.
      *dbl_ptr = S_CAST(double, odds_ratio_numer) / S_CAST(double, odds_ratio_denom);
      return 0;
    }
    // Could use repeated squaring (i.e. port ddr_npwr()) instead if
    // ln_odds_ratio not precomputed.
    dd_real ln_odds_ratio_ddr = *ln_odds_ratio_ddr_ptr;
    if (ln_odds_ratio_ddr.x[0] == DBL_MAX) {
      // todo: error analysis of ddr_sloppy_div.  Not using it for now since
      // this operation does look vulnerable to catastrophic cancellation.
      const dd_real odds_ratio_ddr = ddr_accurate_div(ddr_makei(odds_ratio_numer), ddr_makei(odds_ratio_denom));
      ln_odds_ratio_ddr = ddr_log(odds_ratio_ddr);
      *ln_odds_ratio_ddr_ptr = ln_odds_ratio_ddr;
    }
    *dbl_ptr = exp(ddr_muld(ln_odds_ratio_ddr, odds_ratio_pow).x[0]);
    return 0;
  }
  // possible todo: tune this threshold.
  if (numer_term_ct + denom_term_ct > 256) {
    dd_real ln_odds_ratio_ddr = *ln_odds_ratio_ddr_ptr;
    if (ln_odds_ratio_ddr.x[0] == DBL_MAX) {
      const dd_real odds_ratio_ddr = ddr_accurate_div(ddr_makei(odds_ratio_numer), ddr_makei(odds_ratio_denom));
      ln_odds_ratio_ddr = ddr_log(odds_ratio_ddr);
      *ln_odds_ratio_ddr_ptr = ln_odds_ratio_ddr;
    }
    dd_real starting_lnprobv_ddr = *starting_lnprobv_ddr_ptr;
    if (starting_lnprobv_ddr.x[0] == DBL_MAX) {
      starting_lnprobv_ddr = ddr_muld(ln_odds_ratio_ddr, numer_odds_ratio_pow);
      for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
        starting_lnprobv_ddr = ddr_sub(starting_lnprobv_ddr, ddr_lfact(u31tod(numer_factorial_args[ffac_idx])));
      }
      *starting_lnprobv_ddr_ptr = starting_lnprobv_ddr;
    }
    dd_real lnprobv_ddr = ddr_muld(ln_odds_ratio_ddr, numer_odds_ratio_pow + odds_ratio_pow);
    for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
      lnprobv_ddr = ddr_sub(lnprobv_ddr, ddr_lfact(u31tod(denom_factorial_args[ffac_idx])));
    }
    const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    // ddr_lfact() result has >= 106 bits of precision, should be accurate to
    // 96+ bits (mostly dependent on log).
    // log((2^31)!) is less than 2^36, so we should have at least 60 bits past
    // the decimal point.
    const double epsilon = k2m60 * u31tod(ffac_ct);
    if (lnprob_diff > epsilon) {
      *cmp_resultp = 1;
      *dbl_ptr = exp(lnprob_diff);
      return 0;
    }
    if (lnprob_diff < -epsilon) {
      *cmp_resultp = -1;
      *dbl_ptr = exp(lnprob_diff);
      return 0;
    }
  }
  // numer: Usually CeilPow2(numer_term_ct) + numer_odds_limb_ct limbs in
  //        32-bit limb case.  Also ensure this is at least numer_term_ct + 1,
  //        to cover mpn_mul() and lshift_multilimb() edge cases.
  //        Usually DivUp(CeilPow2(numer_term_ct), 2) + numer_odds_limb_ct
  //        limbs in 64-bit limb case.  Also ensure this is at least
  //        DivUp(numer_term_ct, 2) + 1 to cover edge cases.
  // denom: Similar to above.
  // new_ffac: max(CeilPow2(max_ffac_size), numer_odds_limb_ct,
  //           denom_odds_limb_ct) in 32-bit limb case, DivUp(., 2) in 64-bit
  //           case.
  // generic_wkspace: Used as intermediate buffer for falling-factorial
  //                  calculation, and temporary buffer to copy previous
  //                  numerator or denominator into before multiplication.
  //                  Safe to make this the larger of the numer and denom
  //                  sizes; could tighten this later.

  uint64_t numer_odds_limb_ct = 0;
  uint64_t denom_odds_limb_ct = 0;
  if (odds_ratio_pow && ((odds_ratio_numer != 1) || (odds_ratio_denom != 1))) {
    if (odds_ratio_pow < 0) {
      odds_ratio_pow = -odds_ratio_pow;
      // we don't touch ln_odds_ratio_ddr in this branch, only need to swap
      // odds_ratio_{numer,denom}.
      swap_i64(&odds_ratio_numer, &odds_ratio_denom);
    }
    numer_odds_limb_ct = DivUpU64(bsru64(odds_ratio_numer) * odds_ratio_pow, 32 * kInt32PerLimb);
    denom_odds_limb_ct = DivUpU64(bsru64(odds_ratio_denom) * odds_ratio_pow, 32 * kInt32PerLimb);
  }

  const uint64_t numer_bound1 = numer_odds_limb_ct + (numer_term_ct? DivUp(CeilPow2(numer_term_ct), kInt32PerLimb) : 0);
  uint64_t numer_bound2 = 1 + numer_term_ct;
  const uint64_t denom_bound1 = denom_odds_limb_ct + (denom_term_ct? DivUp(CeilPow2(denom_term_ct), kInt32PerLimb) : 0);
  uint64_t denom_bound2 = 1 + denom_term_ct;
  const uint64_t numer_limb_req = MAXV(numer_bound1, numer_bound2);
  const uint64_t denom_limb_req = MAXV(denom_bound1, denom_bound2);
  const uint64_t new_ffac_limb_req = MAXV(DivUp(CeilPow2(max_ffac_size), kInt32PerLimb), MAXV(numer_odds_limb_ct, denom_odds_limb_ct));
  const uint64_t generic_wkspace_limb_req = MAXV(numer_limb_req, denom_limb_req);
  if (unlikely(generic_wkspace_limb_req > UINT32_MAX)) {
    return 1;
  }
  const uint64_t total_limb_req = numer_limb_req + denom_limb_req + new_ffac_limb_req + generic_wkspace_limb_req;
#ifndef __LP64__
  if (unlikely(total_limb_req > 0x3fffffff)) {
    return 1;
  }
#endif
  if (total_limb_req > *gmp_wkspace_limb_ctp) {
    free_cond(*gmp_wkspacep);
    uintptr_t new_limb_ct = 2 * (*gmp_wkspace_limb_ctp);
    if (new_limb_ct < total_limb_req) {
      new_limb_ct = total_limb_req;
    }
    *gmp_wkspacep = S_CAST(mp_limb_t*, malloc(total_limb_req * sizeof(mp_limb_t)));
    if (unlikely(!(*gmp_wkspacep))) {
      return 1;
    }
    *gmp_wkspace_limb_ctp = new_limb_ct;
  }
  mp_limb_t* numer = *gmp_wkspacep;
  mp_limb_t* denom = &(numer[numer_limb_req]);
  mp_limb_t* new_ffac = &(denom[denom_limb_req]);
  mp_limb_t* generic_wkspace = &(new_ffac[new_ffac_limb_req]);

  uint32_t numer_limb_ct = 1;
  uint32_t denom_limb_ct = 1;
  int64_t pow2 = 0;
  numer[0] = 1;
  denom[0] = 1;
  for (uint32_t ffac_idx = 0; ffac_idx != ffac_ct; ++ffac_idx) {
    const uint32_t numer_factorial_arg = numer_factorial_args[ffac_idx];
    const uint32_t denom_factorial_arg = denom_factorial_args[ffac_idx];
    if (numer_factorial_arg == denom_factorial_arg) {
      continue;
    }
    mp_limb_t* main;
    uint32_t* main_limb_ct_ptr;
    uint32_t top;
    uint32_t ct;
    int32_t sign;
    if (numer_factorial_arg > denom_factorial_arg) {
      main = numer;
      main_limb_ct_ptr = &numer_limb_ct;
      top = numer_factorial_arg;
      ct = top - denom_factorial_arg;
      sign = 1;
    } else {
      main = denom;
      main_limb_ct_ptr = &denom_limb_ct;
      top = denom_factorial_arg;
      ct = top - numer_factorial_arg;
      sign = -1;
    }
    uint32_t new_ffac_limb_ct;
    pow2 += sign * falling_factorial(top, ct, new_ffac, &new_ffac_limb_ct, generic_wkspace);
    const uint32_t main_limb_ct = *main_limb_ct_ptr;
    memcpy(generic_wkspace, main, main_limb_ct * sizeof(mp_limb_t));
    mp_limb_t msl;
    if (new_ffac_limb_ct >= main_limb_ct) {
      msl = mpn_mul(main, new_ffac, new_ffac_limb_ct, generic_wkspace, main_limb_ct);
    } else {
      msl = mpn_mul(main, generic_wkspace, main_limb_ct, new_ffac, new_ffac_limb_ct);
    }
    *main_limb_ct_ptr = main_limb_ct + new_ffac_limb_ct - (msl == 0);
  }
  if (odds_ratio_pow) {
    if (odds_ratio_numer & (odds_ratio_numer - 1)) {
      mul_by_u63_to_u31_pow(odds_ratio_numer, odds_ratio_pow, numer, &numer_limb_ct, new_ffac, generic_wkspace);
    } else {
      // perfect power of 2
      pow2 += odds_ratio_pow * bsru64(odds_ratio_numer);
    }
    if (odds_ratio_denom & (odds_ratio_denom - 1)) {
      mul_by_u63_to_u31_pow(odds_ratio_denom, odds_ratio_pow, denom, &denom_limb_ct, new_ffac, generic_wkspace);
    } else {
      pow2 -= odds_ratio_pow * bsru64(odds_ratio_denom);
    }
  }
  if (pow2 > 0) {
    lshift_multilimb(pow2, numer, &numer_limb_ct);
  } else if (pow2 < 0) {
    lshift_multilimb(-pow2, denom, &denom_limb_ct);
  }
  double numer_d = numer[numer_limb_ct - 1];
  double denom_d = denom[denom_limb_ct - 1];
#ifdef __LP64__
  numer_d *= k2p64;
  if (numer_limb_ct > 1) {
    numer_d += numer[numer_limb_ct - 2];
  }
  denom_d *= k2p64;
  if (denom_limb_ct > 1) {
    denom_d += denom[denom_limb_ct - 2];
  }
#else
  numer_d *= 4294967296.0;
  if (numer_limb_ct > 1) {
    numer_d += numer[numer_limb_ct - 2];
    numer_d *= 4294967296.0;
    if (numer_limb_ct > 2) {
      numer_d += numer[numer_limb_ct - 3];
    }
  }
  denom_d *= 4294967296.0;
  if (denom_limb_ct > 1) {
    denom_d += denom[denom_limb_ct - 2];
    denom_d *= 4294967296.0;
    if (denom_limb_ct > 2) {
      denom_d += denom[denom_limb_ct - 3];
    }
  }
#endif
  double ratio = numer_d / denom_d;
  if (numer_limb_ct == denom_limb_ct) {
    *cmp_resultp = mpn_cmp(numer, denom, numer_limb_ct);
  } else {
    const int32_t limb_diff_ct = numer_limb_ct - denom_limb_ct;
    *cmp_resultp = S_CAST(int32_t, limb_diff_ct);
    ratio = scalbn(ratio, limb_diff_ct * mp_bits_per_limb);
  }
  *dbl_ptr = ratio;
  return 0;
}

// - succ_odds_ratio must be initialized to p / (1-p) reduced to lowest terms,
//   where p is the expected success rate.
//
// - starting_lnprobv_ddr is expected to either be initialized to
//     log(succ_odds_ratio^obs_succ / (obs_succ! (obs_tot - obs_succ)!)),
//   or have x[0] initialized to DBL_MAX to indicate that that calculation
//   hasn't happened.  In the latter case, it may be set to the former value if
//   that is needed in the calculation.
//
// - *cmp_resultp is set to positive value if succ has higher likelihood than
//   obs_succ, 0 if identical likelihood, and negative value if lower
//   likelihood.
//
// - Error is returned iff malloc fails.
BoolErr BinomCompare(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, int32_t succ, dd_real* starting_lnprobv_ddr_ptr, dd_real* ln_odds_ratio_ddr_ptr, intptr_t* cmp_resultp, double* dbl_ptr) {
  // Binomial likelihood is
  //
  //        n!        k        n-k
  //     --------- * p  * (1-p)
  //     k! (n-k)!
  //
  //                        k
  //        n!       [  p  ]         n
  //   = --------- * [ --- ]  * (1-p)
  //     k! (n-k)!   [ 1-p ]
  //
  // where k = # of successes and p is the expected success rate.
  //
  // Thus, the likelihood ratio of interest is
  //
  //   obs_succ! (obs_tot - obs_succ)!                  succ - obs_succ
  //   ------------------------------- * succ_odds_ratio
  //       succ! (obs_tot - succ)!

  uint32_t numer_factorial_args[2];
  numer_factorial_args[0] = obs_succ;
  numer_factorial_args[1] = obs_tot - obs_succ;
  uint32_t denom_factorial_args[2];
  denom_factorial_args[0] = succ;
  denom_factorial_args[1] = obs_tot - succ;

  mp_limb_t* gmp_wkspace = nullptr;
  uintptr_t gmp_wkspace_limb_ct = 0;
  BoolErr reterr = CompareFactorialProductsEx(2, succ_odds_ratio_numer, succ_odds_ratio_denom, succ - obs_succ, obs_succ, numer_factorial_args, denom_factorial_args, starting_lnprobv_ddr_ptr, ln_odds_ratio_ddr_ptr, &gmp_wkspace, &gmp_wkspace_limb_ct, cmp_resultp, dbl_ptr);
  free_cond(gmp_wkspace);
  return reterr;
}

// obs_tot assumed to be <2^31.  succ_odds_ratio_numer and
// succ_odds_ratio_denom must be positive, reduced to lowest terms, have sum <
// 2^63, and represent p/(1-p).
BoolErr BinomLnP(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, uint32_t midp, double* resultp) {
  if (!obs_tot) {
    *resultp = midp? -kLn2 : 0;
    return 0;
  }
  double rate_mult_incr = 1;
  double rate_mult_decr = 1;
  if ((succ_odds_ratio_numer != 1) || (succ_odds_ratio_denom != 1)) {
    // make sure these numbers aren't off by more than 0.5 ULP, even when
    // numerator and/or denominator > 2^53.
    const dd_real numer_ddr = ddr_makei(succ_odds_ratio_numer);
    const dd_real denom_ddr = ddr_makei(succ_odds_ratio_denom);
    rate_mult_incr = ddr_accurate_div(numer_ddr, denom_ddr).x[0];
    rate_mult_decr = ddr_accurate_div(denom_ddr, numer_ddr).x[0];
  }
  double succ = obs_succ;
  double fail = obs_tot - obs_succ;
  // Normalize: succ <= mode.
  if (succ > fail * rate_mult_incr) {
    obs_succ = obs_tot - obs_succ;
    swap_f64(&succ, &fail);
    swap_f64(&rate_mult_incr, &rate_mult_decr);
    swap_i64(&succ_odds_ratio_numer, &succ_odds_ratio_denom);
  }
  const double first_inward_mult = fail * rate_mult_incr / (succ + 1);
  if (!midp) {
    // Might we be at the mode?
    if (first_inward_mult <= 1 + 2 * k2m52) {
      if (first_inward_mult <= 1 - 2 * k2m52) {
        *resultp = 0;
        return 0;
      }
      uint64_t numer_hi;
      uint64_t numer_lo = multiply64to128(obs_tot - obs_succ, succ_odds_ratio_numer, &numer_hi);
      uint64_t denom_hi;
      uint64_t denom_lo = multiply64to128(obs_succ + 1, succ_odds_ratio_denom, &denom_hi);
      if ((denom_hi > numer_hi) || ((denom_hi == numer_hi) && (denom_lo >= numer_lo))) {
        *resultp = 0;
        return 0;
      }
    }
  }
  double lastp = 1;
  double tailp = 1;
  // Iterate outward to floating-point precision limit.
  while (1) {
    fail += 1;
    lastp *= rate_mult_decr * succ / fail;
    succ -= 1;
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
  }
  // In the common case, where we're close enough to the mode that float64
  // underflow/overflow isn't an issue, use the original algorithm: sum all
  // center relative-likelihoods, sum far-tail relative-likelihoods to
  // floating-point precision limit, return log(tailp / (tailp + centerp)).
  //
  // Unfortunately, an extremal rate (e.g. (2^63 - 2)/(2^63 - 1), the maximum
  // possible permitted value) makes center-sum overflow possible much closer
  // to the mode than is the case for the Fisher/HWE exact tests; 17 steps
  // could be enough.  So instead of checking whether we're a constant number
  // of steps from the mode, we compute a lower bound on the number of
  // non-overflowing inward steps we can take from the log of the
  // first-inward-step multiplier (subsequent steps have smaller multipliers).
  succ = obs_succ;
  fail = obs_tot - obs_succ;
  const double ln_mult = log(first_inward_mult);
  double overflow_steps_lower_bound = 0x7fffffff;
  // log(DBL_MAX / 0x7fffffff) = 688.295...
  if (ln_mult > (688.295 / S_CAST(double, 0x7fffffff))) {
    overflow_steps_lower_bound = 688.295 / ln_mult;
  }
  int32_t tie_ct = 1;
  if (succ + overflow_steps_lower_bound > (fail - overflow_steps_lower_bound) * rate_mult_incr) {
    dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
    dd_real ln_odds_ratio_ddr = {{DBL_MAX, 0.0}};
    double one_plus_scaled_eps = 1;
    double centerp = 0;
    lastp = 1;
    while (1) {
      succ += 1;
      lastp *= rate_mult_incr * fail / succ;
      fail -= 1;
      // rate_mult_incr is off by up to 0.5 ULP, and we have two multiplies and
      // a divide.
      one_plus_scaled_eps += 2 * k2m52;
      if (lastp < one_plus_scaled_eps) {
        if (lastp <= 2 - one_plus_scaled_eps) {
          tailp += lastp;
          break;
        }
        // Near-tie.  True value of lastp can be greater than, equal to, or
        // less than 1.
        intptr_t cmp_result;
        if (unlikely(BinomCompare(obs_succ, obs_tot, succ_odds_ratio_numer, succ_odds_ratio_denom, S_CAST(int32_t, succ), &starting_lnprobv_ddr, &ln_odds_ratio_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        one_plus_scaled_eps = 1 + 3 * k2m52;
        if (cmp_result <= 0) {
          tailp += lastp;
          tie_ct += (cmp_result == 0);
          break;
        }
      }
      centerp += lastp;
    }
    // TODO
  }
  // TODO
  return 0;
}

#ifdef __cplusplus
}
#endif
