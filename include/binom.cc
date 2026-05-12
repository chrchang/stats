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

// Assumes n < 2^52.
double LnBinomCoeff(int64_t n, int64_t k) {
  if ((k == 0) || (k == n)) {
    return 0;
  }
  return ddr_sub(ddr_lfact(n),
                 ddr_add_lfacts(k, n-k)).x[0];
}

// Currently assumes k < n.
dd_real binom_ln_prob_internal(int64_t k, int64_t n, dd_real p_ddr) {
  const dd_real q_ddr = ddr_sub(ddr_maked(1), p_ddr);
  dd_real ln_q_ddr = ddr_negate(_ddr_log2);
  if ((p_ddr.x[0] != 0.5) || (p_ddr.x[1] != 0.0)) {
    ln_q_ddr = ddr_log(q_ddr);
  }
  const dd_real nmk_ln_q_ddr = ddr_muld(ln_q_ddr, n-k);
  if (k == 0) {
    return nmk_ln_q_ddr;
  }
  dd_real ln_p_ddr = ddr_negate(_ddr_log2);
  if ((p_ddr.x[0] != 0.5) || (p_ddr.x[1] != 0.0)) {
    ln_p_ddr = ddr_log(p_ddr);
  }
  const dd_real k_ln_p_ddr = ddr_muld(ln_p_ddr, k);
  dd_real ddrs[5];
  ddrs[0] = k_ln_p_ddr;
  ddrs[1] = nmk_ln_q_ddr;
  ddrs[2] = ddr_lfact(n);
  ddrs[3] = ddr_negate(ddr_lfact(k));
  ddrs[4] = ddr_negate(ddr_lfact(n-k));
  return ddr_sort_and_add(5, ddrs);
}

// Assumes 0 <= k <= n < 2^52, 0 < p < 1.
// If p is too close to 1 to be well-represented by a dd_real, pass in (n-k, n,
// 1-p) instead.
double BinomMass(int64_t k, int64_t n, dd_real p_ddr, uint32_t logp) {
  if (k == n) {
    k = 0;
    p_ddr = ddr_sub(ddr_maked(1), p_ddr);
  }
  const dd_real ln_prob_ddr = binom_ln_prob_internal(k, n, p_ddr);
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  // Note that if we want to deliver j-bit precision for x < 2^{-512}, log(x)
  // needs to be represented to ~(j+9)-bit precision.
  return ddr_exp(ln_prob_ddr).x[0];
}

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
// - ln_odds_ratio_ddr is expected to either be initialized to log(odds_ratio),
//   or have x[0] initialized to DBL_MAX, etc.
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
// succ_odds_ratio_{numer,denom} must be positive, reduced to lowest terms,
// have sum < 2^63, and represent p/(1-p).
// Only the main loops are speed-optimized for now.  There is some setup
// overhead for odds_ratio != 1 which is currently written to avoid
// proliferation of special cases, but can be easily accelerated if we know
// e.g. numerator and denominator < 2^31.
BoolErr BinomLnP(int32_t obs_succ, int32_t obs_tot, int64_t succ_odds_ratio_numer, int64_t succ_odds_ratio_denom, int32_t midp, double* resultp) {
  if (!obs_tot) {
    *resultp = midp? -kLn2 : 0;
    return 0;
  }
  double succ = obs_succ;
  double fail = obs_tot - obs_succ;
  // Normalize: succ <= mode.
  // (even if there is rounding error, this is enough to guarantee that succ-1
  // has lower likelihood than succ.)
  if (succ * S_CAST(double, succ_odds_ratio_denom) > fail * S_CAST(double, succ_odds_ratio_numer)) {
    obs_succ = obs_tot - obs_succ;
    swap_f64(&succ, &fail);
    swap_i64(&succ_odds_ratio_numer, &succ_odds_ratio_denom);
  }
  dd_real succ_odds_ratio_ddr = ddr_maked(1);
  if ((succ_odds_ratio_numer != 1) || (succ_odds_ratio_denom != 1)) {
    // make sure these numbers aren't off by more than 0.5 ULP, even when
    // numerator and/or denominator > 2^53.
    const dd_real numer_ddr = ddr_makei(succ_odds_ratio_numer);
    const dd_real denom_ddr = ddr_makei(succ_odds_ratio_denom);
    succ_odds_ratio_ddr = ddr_accurate_div(numer_ddr, denom_ddr);
  }
  const double succ_odds_ratio = succ_odds_ratio_ddr.x[0];
  const double first_inward_mult = fail * succ_odds_ratio / (succ + 1);
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
  double tailp = 1 - midp * 0.5;
  // Iterate outward to floating-point precision limit.
  while (1) {
    fail += 1;
    lastp *= succ / (succ_odds_ratio * fail);
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
  // non-overflowing inward steps we can take, using the log of the
  // first-inward-step multiplier (subsequent steps have smaller multipliers).
  const double smallest_evaluated_succ = succ;
  succ = obs_succ;
  fail = obs_tot - obs_succ;
  const double ln_mult = log(first_inward_mult);
  double overflow_steps_lower_bound = 0x7fffffff;
  // log(DBL_MAX / (2^31 - 1)) = 688.295...
  if (ln_mult > (688.295 / S_CAST(double, 0x7fffffff))) {
    overflow_steps_lower_bound = 688.295 / ln_mult;
  }
  // succ_odds_ratio * (tot - modal_succ) / modal_succ = 1
  // succ_odds_ratio * (tot - modal_succ) = modal_succ
  // succ_odds_ratio * tot = modal_succ * (1 + succ_odds_ratio)
  // possible for modal_succ to round up to just obs_tot
  const double obs_totd = obs_tot;
  const double modal_succ = obs_totd * succ_odds_ratio / (1 + succ_odds_ratio);
  if (succ + overflow_steps_lower_bound > modal_succ) {
    dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
    dd_real ln_odds_ratio_ddr = {{DBL_MAX, 0.0}};
    double one_plus_scaled_eps = 1 + k2m52;
    double centerp = midp * 0.5;
    lastp = 1;
    while (1) {
      succ += 1;
      lastp *= succ_odds_ratio * fail / succ;
      fail -= 1;
      // succ_odds_ratio is off by up to 0.5 ULP, and we have two multiplies
      // and a divide.
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
          if (midp && (cmp_result == 0)) {
            tailp -= 0.5;
            centerp += 0.5;
          }
          break;
        }
      }
      centerp += lastp;
    }
    // Continue down tail to floating-point precision limit.
    while (1) {
      succ += 1;
      lastp *= succ_odds_ratio * fail / succ;
      fail -= 1;
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
    }
    *resultp = log(tailp / (tailp + centerp));
    return 0;
  }
  dd_real ln_odds_ratio_ddr = ddr_log(succ_odds_ratio_ddr);
  dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(ln_odds_ratio_ddr, succ),
            ddr_add_lfacts(succ, fail));
  dd_real lnfail_ddr = {{-_ddr_log2.x[0], -_ddr_log2.x[1]}};
  if (!ddr_is_zero(ln_odds_ratio_ddr)) {
    // probable todo: this is a bit redundant with earlier initialization
    const dd_real fail_ddr = ddr_makei(succ_odds_ratio_denom);
    const dd_real succ_plus_fail_ddr = ddr_makeu64(S_CAST(uint64_t, succ_odds_ratio_numer) + S_CAST(uint64_t, succ_odds_ratio_denom));
    lnfail_ddr = ddr_log(ddr_accurate_div(fail_ddr, succ_plus_fail_ddr));
  }
  const dd_real lnprobf_ddr =
    ddr_add(ddr_lfact(obs_totd), ddr_muld(lnfail_ddr, obs_totd));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];

  // Now we want to jump near the other tail, without evaluating that many
  // contingency table log-likelihoods along the way.
  //
  // Each full log-likelihood evaluation requires 2 ddr_lfact() calls.  Since
  // they are now performed with extra precision, they require hundreds of
  // floating-point operations, so we want to limit ourselves to 1-2 full
  // evaluations most of the time.  (Possible todo: use lower-accuracy Lfact()
  // to jump around, followed by ddr_lfact() when exiting the loop.  Should be
  // an easy performance win, but there's a complexity cost so I'll wait until
  // I see a scenario where this branch executes frequently...)
  //
  // The current heuristic starts by reflecting (smallest_evaluated_succ +
  // succ) * 0.5 across the (continuous) mode, performing a full log-likelihood
  // check at the nearest valid point.  Hopefully we find that we're in
  // (starting_lnprob - tolerance, starting_lnprob], so we're at or near a
  // table that actually contributes to the tail-sum; unlike the Fisher's and
  // HWE cases, we can't fix the tolerance at 62 * kLn2, but we can compute a
  // value >= 53 * kLn2 large enough to guarantee at least one point falls
  // inside.
  //
  // If not, we jump again, using Newton's method.
  // If succ is too low (i.e. current log-likelihood is too high), when we
  // increase succ by 1, the likelihood gets multiplied by
  //   succ_odds_ratio * fail / (succ+1)
  // i.e. we're adding the logarithm of this value to the log-likelihood.
  // If succ is too high, when we decrease succ by 1, the likelihood gets
  // multiplied by
  //   succ / (succ_odds_ratio * (fail+1))
  // We use the log of the first expression as the Newton's method f'(x) when
  // we're jumping to higher succ, and the negative-log of the second
  // expression when we're jumping to lower homr.
  // f''(x) is always negative, so we can aim for starting_lnprob instead of
  // the middle of the interval.

  // L(obs_tot) / L(obs_tot-1) = succ_odds_ratio * 1 / obs_tot
  // If this value is >= 1, obs_tot is a mode.  Separating out that case makes
  // the remaining logic simpler.
  if (ddr_geqd(succ_odds_ratio_ddr, obs_totd)) {
    *resultp = starting_lnprob + log(tailp);
    return 0;
  }

  succ = 2 * modal_succ - (succ + smallest_evaluated_succ) * 0.5;
  if (succ > obs_totd) {
    succ = obs_totd;
  }
  succ = S_CAST(int32_t, succ);

  // obs_tot is past the mode, and |log(L(obs_tot) / L(obs_tot-1))| is the
  // largest gap between adjacent log-likelihoods on this tail.  Set
  // |lnprob_diff_min| >= this value.
  double lnprob_diff_min = log(succ_odds_ratio / obs_totd) * (1 + kSmallEpsilon);
  if (lnprob_diff_min > -53 * kLn2) {
    lnprob_diff_min = -53 * kLn2;
  }

  while (1) {
    fail = obs_totd - succ;
    const dd_real lnprobv_ddr =
      ddr_sub(ddr_muld(ln_odds_ratio_ddr, succ),
              ddr_add_lfacts(succ, fail));
    const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    if (lnprob_diff >= k2m53) {
      if (fail == 0) {
        *resultp = starting_lnprob + log(tailp);
        return 0;
      }
      const double ll_deriv = ln_odds_ratio_ddr.x[0] + log(fail / (succ + 1));
      // Round up, to guarantee that we make progress.
      // (lnprob_diff is positive and ll_deriv is negative.)
      // This may overshoot.  But the function is guaranteed to terminate
      // because we never overshoot (and we do always make progress on each
      // step) once we're on the other side.
      succ += ceil_smalleps(-lnprob_diff / ll_deriv);
      if (succ > obs_totd) {
        succ = obs_totd;
      }
    } else if (lnprob_diff > lnprob_diff_min) {
      lastp = exp(lnprob_diff);
      break;
    } else {
      const double ll_deriv = ln_odds_ratio_ddr.x[0] + log((fail + 1) / succ);
      // Round down, to guarantee we don't overshoot.
      // |lnprob_diff| >= |lnprob_diff_min| > |ll_deriv| so we're guaranteed
      // to make progress.
      succ -= S_CAST(int64_t, lnprob_diff / ll_deriv);
    }
  }
  // Sum toward center, until lastp >= 1.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double lastp_tail = lastp;
  const double succ_tail = succ;
  while (lastp <= one_minus_scaled_eps) {
    tailp += lastp;
    fail += 1;
    lastp *= succ / (succ_odds_ratio * fail);
    succ -= 1;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lastp < 2 - one_minus_scaled_eps) {
    intptr_t cmp_result;
    if (unlikely(BinomCompare(obs_succ, obs_tot, succ_odds_ratio_numer, succ_odds_ratio_denom, S_CAST(int32_t, succ), &starting_lnprobv_ddr, &ln_odds_ratio_ddr, &cmp_result, &lastp))) {
      return 1;
    }
    if (cmp_result <= 0) {
      tailp += lastp;
      if (midp && (cmp_result == 0)) {
        tailp -= 0.5;
      }
    }
  }
  // Sum away from center, until sums stop changing.
  lastp = lastp_tail;
  succ = succ_tail;
  fail = obs_totd - succ;
  while (1) {
    succ += 1;
    lastp *= succ_odds_ratio * fail / succ;
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
    fail -= 1;
  }
  *resultp = starting_lnprob + log(tailp);
  return 0;
}

// ibeta_fraction2_ln_eps() and dependencies adapted from Boost 1.91.0.  This
// derived code is subject to the following license:
//
// *****
// Boost Software License - Version 1.0 - August 17th, 2003
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
// *****

static const double kLentzFpmin = DBL_MIN * 16;

// We want more than 6 digits of accuracy here.
static const double kLanczosDoubleSumDenom[13] = {0, 39916800, 120543840, 150917976, 105258076, 45995730, 13339535, 2637558, 357423, 32670, 1925, 66, 1};
static const double kLanczosDoubleSumExpgNumer[13] = {
  56906521.91347156388090791033559122686859,
  103794043.1163445451906271053616070238554,
  86363131.28813859145546927288977868422342,
  43338889.32467613834773723740590533316085,
  14605578.08768506808414169982791359218571,
  3481712.15498064590882071018964774556468,
  601859.6171681098786670226533699352302507,
  75999.29304014542649875303443598909137092,
  6955.999602515376140356310115515198987526,
  449.9445569063168119446858607650988409623,
  19.51992788247617482847860966235652136208,
  0.5098416655656676188125178644804694509993,
  0.006061842346248906525783753964555936883222
};

// this depends on the polynomial coefficients above
// exactly 808618867 * 2^{-27}, don't need to represent this as dd_real
static const double kLanczosDoubleG = 6.024680040776729583740234375;

double lanczos_sum_d_expg_scaled_imp(double zz, double* s2_ptr) {
  double s1;
  double s2;
  // zz currently guaranteed to be >1.
  /*
  if (zz <= 1) {
    s1 = kLanczosDoubleSumExpgNumer[12];
    s2 = kLanczosDoubleSumDenom[12];
    for (int32_t ii = 11; ii >= 0; --ii) {
      s1 *= zz;
      s2 *= zz;
      s1 += kLanczosDoubleSumExpgNumer[S_CAST(uint32_t, ii)];
      s2 += kLanczosDoubleSumDenom[S_CAST(uint32_t, ii)];
    }
  } else {
  */
  zz = 1 / zz;
  s1 = kLanczosDoubleSumExpgNumer[0];
  s2 = kLanczosDoubleSumDenom[0];
  for (uint32_t uii = 1; uii != 13; ++uii) {
    s1 *= zz;
    s2 *= zz;
    s1 += kLanczosDoubleSumExpgNumer[uii];
    s2 += kLanczosDoubleSumDenom[uii];
  }
  // }
  *s2_ptr = s2;
  return s1;
}

dd_real ibeta_power_terms_d_ln(double aa, double bb, dd_real p_ddr, dd_real q_ddr, dd_real ay_minus_bx_ddr) {
  // returns log((x^a)(y^b) / Beta(a,b))
  //
  // normalized always true
  // prefix always 1
  // aa and bb always large
  double cc = aa + bb;
  const dd_real gh_ddr = {{kLanczosDoubleG - 0.5, 0.0}};
  const dd_real agh_ddr = ddr_addd(gh_ddr, aa);
  const dd_real bgh_ddr = ddr_addd(gh_ddr, bb);
  const dd_real cgh_ddr = ddr_addd(gh_ddr, cc);

  double numer_a;
  const double denom_a = lanczos_sum_d_expg_scaled_imp(aa, &numer_a);
  double numer_b;
  const double denom_b = lanczos_sum_d_expg_scaled_imp(bb, &numer_b);
  double denom_c;
  const double numer_c = lanczos_sum_d_expg_scaled_imp(cc, &denom_c);
  // Calculate result with ordinary precision; pointless to go further here
  // unless we widen Lanczos sum calculations.
  double result = (numer_a * numer_b * numer_c) / (denom_a * denom_b * denom_c);
  result *= sqrt(agh_ddr.x[0] * bgh_ddr.x[0] * kRecipE / cgh_ddr.x[0]);
  // Represent log-probability with more than 53 bits, to preserve non-log
  // accuracy without sacrificing range.
  // TODO: Can we prove <this result> / <ibeta_fraction2's ff> doesn't
  // overflow or underflow?  Then we can replace one of these expensive ddr_log
  // operations with an ordinary multiply.
  // todo: ddr_expd, ddr_logd?
  dd_real result_ln_ddr = ddr_log(ddr_maked(result));
  // Calculate l1 and l2 with extra precision, since magnitude can greatly
  // exceed that of ln(result).
  const dd_real l1_ddr = ddr_accurate_div(ddr_negate(ddr_add(ay_minus_bx_ddr, ddr_muld(q_ddr, gh_ddr.x[0]))), agh_ddr);
  const dd_real l2_ddr = ddr_accurate_div(ddr_sub(ay_minus_bx_ddr, ddr_muld(p_ddr, gh_ddr.x[0])), bgh_ddr);
  return ddr_add3(result_ln_ddr,
                  ddr_muld(ddr_log1p(l1_ddr), aa),
                  ddr_muld(ddr_log1p(l2_ddr), bb));
}

dd_real ibeta_fraction2_ln_ddr(double aa, double bb, dd_real p_ddr, dd_real q_ddr, uint32_t inv) {
  // normalized always true, min(aa, bb) >= 40, max much larger

  const dd_real ay_minus_bx_ddr = ddr_sub(ddr_muld(q_ddr, aa), ddr_muld(p_ddr, bb));
  dd_real result_ln_ddr = ibeta_power_terms_d_ln(aa, bb, p_ddr, q_ddr, ay_minus_bx_ddr);

  // see Boost continued_fraction_b()
  const double ay_minus_bx_plus1 = ay_minus_bx_ddr.x[0] + 1.0;
  double ff = (aa * ay_minus_bx_plus1) / (aa + 1.0);
  if (ff == 0.0) {
    ff = kLentzFpmin;
  }
  const double xx = p_ddr.x[0];
  const double x2 = xx * xx;
  const double two_minus_x = 2 - xx;
  double cc = ff;
  double dd = 0.0;
  double mm = 1.0;
  while (1) {
    const double denom = aa + 2 * mm - 1;
    const double cur_a = (mm * (aa + mm - 1) / denom) * ((aa + bb + mm - 1) / denom) * (bb - mm) * x2;
    double cur_b = mm;
    cur_b += (mm * (bb - mm) * xx) / denom;
    cur_b += ((aa + mm) * (ay_minus_bx_plus1 + mm * two_minus_x)) / (aa + 2 * mm + 1);
    mm += 1.0;
    dd = cur_b + cur_a * dd;
    if (dd == 0.0) {
      dd = kLentzFpmin;
    }
    cc = cur_b + cur_a / cc;
    if (cc == 0.0) {
      cc = kLentzFpmin;
    }
    dd = 1.0 / dd;
    const double delta = cc * dd;
    if (delta == 1.0) {
      result_ln_ddr = ddr_addd(result_ln_ddr, -log(ff));
      if (!inv) {
        return result_ln_ddr;
      }
      return ddr_log1p(ddr_negate(ddr_exp(result_ln_ddr)));
    }
    ff *= delta;
  }
}

// Requires 0 <= obs_succ <= obs_tot < 2^52 and 2^{-960} < p < 1.
// (Sometimes works for 0 < p <= 2^{-960}, but let's leave that out of the
// function contract until the holes in that region are plugged in.)
//
// Larger obs_succ and obs_tot are allowed than for the 2-sided test because
// there's no risk of needing to expand gigantic ratios of factorials to handle
// likelihood near-ties correctly.
//
// Benchmark results revealed that Boost 1.91 ibetac(k+1, n-k, p) (which is
// called by scipy.stats.binom.logcdf()) became faster than this function's
// initial implementation once obs_tot was ~1000, and its results were
// acceptably accurate.  So we now use its main algorithm when obs_tot > 2^10
// and min(k+1, n-k) >= 40.
//
// Interestingly, the scipy implementation has much higher overhead, even after
// initialization, despite relying on nearly identical C++ code.  E.g.
//
//   >>> import exact_tests, scipy, timeit
//   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
//   0.017715082969516516
//   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
//   0.01739558414556086
//   >>> timeit.timeit(lambda: exact_tests.pbinom(157000000, 419430500, 0.375, approx=True), number=10000)
//   0.018672874895855784
//   >>> timeit.timeit(lambda: scipy.stats.binom.logcdf(157000000, 419430500, 0.375), number=10000)
//   1.005605333019048
//   >>> timeit.timeit(lambda: scipy.stats.binom.logcdf(157000000, 419430500, 0.375), number=10000)
//   0.22769695916213095
//   >>> timeit.timeit(lambda: scipy.stats.binom.logcdf(157000000, 419430500, 0.375), number=10000)
//   0.24197920807637274
//   >>> timeit.timeit(lambda: scipy.stats.binom.cdf(157000000, 419430500, 0.375), number=10000)
//   0.27701791608706117
//   >>> timeit.timeit(lambda: scipy.stats.binom.cdf(157000000, 419430500, 0.375), number=10000)
//   0.27656474988907576
//   >>> timeit.timeit(lambda: scipy.stats.binom.cdf(157000000, 419430500, 0.375), number=10000)
//   0.27791083394549787

// See Pbinom() below for a higher-accuracy variant of this function; this one
// is limited by the float64 precision of the Lanczos sums, and everything else
// here is tuned to that level of relative error.
double BinomOneSidedLnP(int64_t obs_succ, int64_t obs_tot, dd_real p_ddr, uint32_t succ_is_greater_alt, int32_t midp, uint32_t logp) {
  // Normalize to alternative hypothesis = #succ-less-than-expected, so the
  // remainder is essentially a cmf calculation with parameter succ.
  dd_real q_ddr = ddr_negate(ddr_addd(p_ddr, -1.0));
  if (succ_is_greater_alt) {
    obs_succ = obs_tot - obs_succ;
    swap_ddr(&p_ddr, &q_ddr);
  }
  if ((obs_tot > 1024) && (MINV(obs_succ, obs_tot - obs_succ) >= 40)) {
    double aa = obs_succ + 1;
    double bb = obs_tot - obs_succ;
    uint32_t inv = 1;
    double lambda;
    if (aa < bb) {
      lambda = prefer_fma(aa + bb, -p_ddr.x[0], aa);
    } else {
      lambda = prefer_fma(aa + bb, q_ddr.x[0], -bb);
    }
    if (lambda < 0.0) {
      swap_f64(&aa, &bb);
      swap_ddr(&p_ddr, &q_ddr);
      inv = !inv;
    }
    // TODO: use_asym branch for gigantic cases
    /*
    uint32_t use_asym = 0;
    const double ma = MAXV(aa, bb);
    const double xa = (ma == aa)? xx : yy;
    const double saddle = ma / (aa + bb);
    if ((ma > (0.00001 / k2m53)) && (ma / MINV(aa, bb) < ((xa < saddle)? 2 : 15))) {
      if (aa == bb) {
        use_asym = 1;
      } else {
        double powers = exp(log(xx / (aa / (aa + bb))) * aa + log(yy / (bb / (aa + bb))) * bb);
        if (powers < k2m53) {
          use_asym = 1;
        }
      }
    }
    */

    dd_real result_ln_ddr = ibeta_fraction2_ln_ddr(aa, bb, p_ddr, q_ddr, inv);
    if (midp) {
      // Subtract 0.5 * pmf(k, n, p).
      const dd_real ln_half_pmf_ddr = ddr_sub(binom_ln_prob_internal(obs_succ, obs_tot, p_ddr), _ddr_log2);
      const dd_real ln_ratio_ddr = ddr_sub(ln_half_pmf_ddr, result_ln_ddr);
      result_ln_ddr = ddr_add(result_ln_ddr, ddr_log(ddr_negate(ddr_expm1(ln_ratio_ddr))));
    }
    if (logp) {
      return result_ln_ddr.x[0];
    }
    return ddr_exp(result_ln_ddr).x[0];
  }
  const double succ_odds_ratio = ddr_accurate_div(p_ddr, q_ddr).x[0];
  double succ = obs_succ;
  double fail = obs_tot - obs_succ;
  if (succ > fail * succ_odds_ratio) {
    // We're at or to the right of the mode.
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know the p-value > 1 - DBL_MIN (at which point
    // we just return 0, don't want to risk imposing a surprising
    // denormal-handling performance penalty for no good reason), or remaining
    // left likelihoods are smaller than the precision limit.
    const double first_right_mult = succ_odds_ratio * fail / (succ + 1);
    // r + r^2 + ... = r / (1-r)
    const double right_upper_bound = 0.5 * midp + first_right_mult / (1 - first_right_mult);
    if (right_upper_bound == 0.0) {
      // p-value is exactly 1 when fail==0 and midp is false
      return logp? 0 : 1;
    }

    // Scale our starting likelihood so that we overflow to INFINITY when we'd
    // want to early-exit and return 0; this saves us a comparison in the loop.
    const double startp = (DBL_MAX * (logp? DBL_MIN : (1.0 / (1LL << 54)))) / right_upper_bound;
    double lastp = startp;
    double left_sum = startp;
    while (1) {
      fail += 1;
      lastp *= succ / (succ_odds_ratio * fail);
      succ -= 1;
      const double preaddp = left_sum;
      left_sum += lastp;
      if (left_sum == preaddp) {
        break;
      }
    }
    if (left_sum == INFINITY) {
      return logp? 0 : 1;
    }
    // startp is now potentially < 1, so this operation might throw away some
    // range.
    // left_sum /= startp;

    // Now compute the right-sum to the precision limit.
    double right_sum = first_right_mult * startp;
    succ = obs_succ + 1;
    fail = obs_tot - obs_succ - 1;
    lastp = right_sum;
    while (1) {
      succ += 1;
      lastp *= succ_odds_ratio * fail / succ;
      fail -= 1;
      const double preaddp = right_sum;
      right_sum += lastp;
      if (right_sum == preaddp) {
        break;
      }
    }
    // For one-sided test, slightly more convenient to exclude midp term from
    // left_sum and right_sum since it just cancels out in denom
    const double midp_numer = -0.5 * midp * startp;
    const double denom = right_sum + left_sum;
    if (!logp) {
      return (left_sum + midp_numer) / denom;
    }
    return log1p((midp_numer - right_sum) / denom);
  }
  // We're to the left of the mode, and are responsible for tiny p-values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise... if the problem instance isn't *that* large, we could use
  // Lfact() to compute the starting log-likelihood and eat a big catastrophic
  // cancellation error, but I'll start with just the slow-and-accurate
  // calculation.
  const double ln_first_inward_mult = log(succ_odds_ratio * fail / (succ + 1));
  double overflow_steps_lower_bound = 1LL << 52;
  // log(DBL_MAX / 2^52) = 673.739...
  if (ln_first_inward_mult > 673.739 * k2m52) {
    overflow_steps_lower_bound = 673.739 / ln_first_inward_mult;
  }
  const double obs_totd = obs_tot;
  const double modal_succ = obs_totd * succ_odds_ratio / (1 + succ_odds_ratio);
  if (succ + overflow_steps_lower_bound > modal_succ) {
    double lastp = 1;
    double right_sum = 0;
    while (1) {
      succ += 1;
      lastp *= succ_odds_ratio * fail / succ;
      fail -= 1;
      const double preaddp = right_sum;
      right_sum += lastp;
      if (right_sum == preaddp) {
        break;
      }
    }
    succ = obs_succ;
    fail = obs_tot - obs_succ;
    lastp = 1;
    double left_sum = 1;
    while (1) {
      fail += 1;
      lastp *= succ / (succ_odds_ratio * fail);
      succ -= 1;
      const double preaddp = left_sum;
      left_sum += lastp;
      if (left_sum == preaddp) {
        break;
      }
    }
    const double pval = (left_sum - 0.5 * midp) / (left_sum + right_sum);
    return logp? log(pval) : pval;
  }
  const dd_real succ_odds_ratio_ddr = ddr_maked(succ_odds_ratio);
  const dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(ddr_log(succ_odds_ratio_ddr), succ),
            ddr_add_lfacts(succ, fail));
  dd_real lnfail_ddr = {{-_ddr_log2.x[0], -_ddr_log2.x[1]}};
  if (succ_odds_ratio != 1.0) {
    // log(1 / (1 + succ_odds_ratio)) = -log(1 + succ_odds_ratio)
    lnfail_ddr = ddr_negate(ddr_log(ddr_addd(succ_odds_ratio_ddr, 1.0)));
  }
  const dd_real lnprobf_ddr =
    ddr_add(ddr_lfact(obs_totd), ddr_muld(lnfail_ddr, obs_totd));
  const dd_real starting_lnprob_ddr = ddr_add(lnprobf_ddr, starting_lnprobv_ddr);
  double lastp = 1;
  double left_sum = 1 - 0.5 * midp;
  while (1) {
    fail += 1;
    lastp *= succ / (succ_odds_ratio * fail);
    succ -= 1;
    const double preaddp = left_sum;
    left_sum += lastp;
    if (left_sum == preaddp) {
      break;
    }
  }
  if (!logp) {
    // Add log(2^52) to starting_lnprob before exponentiating, so it doesn't
    // underflow unless the final result would also underflow.
    // todo: precompute this
    const dd_real _ddr_52_log2 = ddr_muld(_ddr_log2, 52);
    return left_sum * k2m52 * ddr_exp(ddr_add(starting_lnprob_ddr, _ddr_52_log2)).x[0];
  }
  return starting_lnprob_ddr.x[0] + log(left_sum);
}

// Assumes 0 <= k < n < 2^52, 2^{-960} <= p < 1.
// Should achieve <0.6 ULP relative error except when n is well over 2^31.
//
// See BinomOneSidedLnP() above for the faster variant of this function which
// doesn't try to get the last few bits right.
//
// TODO: What can we prove about a variant of ibeta_fraction2_ln_ddr which
// computes 17- (or 24-?)term Lanczos sums with dd_real precision?  Can we
// safely hand off most huge-n cases to it?
double Pbinom(int64_t obs_k, int64_t n, dd_real p_ddr, uint32_t logp) {
  const dd_real q_ddr = ddr_ieee_sub(ddr_maked(1.0), p_ddr);
  const dd_real pdq_ddr = ddr_accurate_div(p_ddr, q_ddr);
  const dd_real qdp_ddr = ddr_accurate_div(q_ddr, p_ddr);
  double k = obs_k;
  double nmk = n - obs_k;
  if (k > nmk * pdq_ddr.x[0]) {
    // We're at or to the right of the mode.
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know we can safely return a cmf value of 1, or
    // remaining left likelihoods are smaller than the precision limit.
    const dd_real first_right_mult_ddr = ddr_divd(ddr_muld(pdq_ddr, nmk), k+1);
    const double right_upper_bound = first_right_mult_ddr.x[0] / ddr_negate(ddr_subd(first_right_mult_ddr, 1.0)).x[0];
    // Function contract ensures right_upper_bound > 0.

    // Scale our starting likelihood so that we overflow to INFINITY/nan when
    // we'd want to early-exit and return 1; this saves us a comparison in the
    // loop.
    const double startp = (DBL_MAX * (logp? DBL_MIN : (k2m53 * 0.5))) / right_upper_bound;
    // We want to compute left_sum_ddr to at most 2^{-57} relative error, but
    // we also want to drop down from this slow dd_real-based loop to the much
    // faster float64-based loop as soon as we can prove that won't make us
    // miss the accuracy target.
    //
    // Initial lastp error when we exit the first loop will be < 0.501 ULP; the
    // lastp values we add to left_tail_sum have error bounded above by 2.5
    // ULPs, 4.5 ULPs, 6.5 ULPs, etc.  Thus, left_tail_sum should have error
    // bounded above by 3 ULPs, then 8, 15, 24, 35, ...; and k^2 is a loose
    // upper bound on the final error.
    //
    // 2^{-57} corresponds to at least 0.03125 ULPs.  If k is so large that
    // min_mult ends up just being 1, that's ok.
    const double min_mult_left = 1 + 0.03125 / (k * k);
    dd_real lastp_ddr = ddr_maked(startp);
    dd_real left_sum_ddr = lastp_ddr;
    while (1) {
      nmk += 1;
      lastp_ddr = ddr_mul(lastp_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
      k -= 1;
      const double preaddp = left_sum_ddr.x[0];
      left_sum_ddr = ddr_add(left_sum_ddr, lastp_ddr);
      // This overflows to nan rather than INFINITY in my testing; I'll write
      // this to be agnostic to that detail.
      if (!(left_sum_ddr.x[0] > preaddp * min_mult_left)) {
        break;
      }
      if (preaddp != preaddp) {
        return 0.0;
      }
    }
    if (!(left_sum_ddr.x[0] < INFINITY)) {
      return logp? 0.0 : 1.0;
    }
    if (k > 0) {
      // Continue the calculation with ordinary precision.
      const double pdq = pdq_ddr.x[0];
      double lastp = lastp_ddr.x[0];
      double left_tail_sum = 0.0;
      while (1) {
        nmk += 1;
        lastp *= k / (pdq * nmk);
        k -= 1;
        const double preaddp = left_tail_sum;
        left_tail_sum += lastp;
        if (left_tail_sum == preaddp) {
          break;
        }
      }
      left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
      if (left_sum_ddr.x[0] == INFINITY) {
        return logp? 0.0 : 1.0;
      }
    }

    // Now compute the right-sum to at most 2^{-57} relative error.
    dd_real right_sum_ddr = ddr_muld(first_right_mult_ddr, startp);
    k = obs_k + 1;
    nmk = n - obs_k - 1;
    lastp_ddr = right_sum_ddr;
    if (nmk > 0) {
      const double min_mult_right = 1 + 0.03125 / (nmk * nmk);
      while (1) {
        k += 1;
        lastp_ddr = ddr_mul(lastp_ddr, ddr_divd(ddr_muld(pdq_ddr, nmk), k));
        nmk -= 1;
        const double preaddp = right_sum_ddr.x[0];
        right_sum_ddr = ddr_add(right_sum_ddr, lastp_ddr);
        if (right_sum_ddr.x[0] <= preaddp * min_mult_right) {
          break;
        }
      }
      if (nmk > 0) {
        // Continue the calculation with ordinary precision.
        const double pdq = pdq_ddr.x[0];
        double lastp = lastp_ddr.x[0];
        double right_tail_sum = 0.0;
        while (1) {
          k += 1;
          lastp *= pdq * nmk / k;
          nmk -= 1;
          const double preaddp = right_tail_sum;
          right_tail_sum += lastp;
          if (right_tail_sum == preaddp) {
            break;
          }
        }
        right_sum_ddr = ddr_addd(right_sum_ddr, right_tail_sum);
      }
    }
    const dd_real denom_ddr = ddr_add(left_sum_ddr, right_sum_ddr);
    if (!(denom_ddr.x[0] < INFINITY)) {
      return logp? 0.0 : 1.0;
    }
    const dd_real one_minus_prob_ddr = ddr_accurate_div(right_sum_ddr, denom_ddr);
    if (!logp) {
      return -ddr_subd(one_minus_prob_ddr, 1.0).x[0];
    }
    return ddr_log1p(ddr_negate(one_minus_prob_ddr)).x[0];
  }
  // We're to the left of the mode, and are responsible for tiny cmf values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise, we use ddr_lfact() and friends to compute the starting
  // log-likelihood (this is the only step that could introduce >0.1 ULP error,
  // and then only for huge n), and accumulate the tail-sum from there.
  const double ln_first_inward_mult = log(pdq_ddr.x[0] * nmk / (k + 1));
  double overflow_steps_lower_bound = 1LL << 52;
  // log(DBL_MAX / 2^52) = 673.739...
  if (ln_first_inward_mult > 673.739 * k2m52) {
    overflow_steps_lower_bound = 673.739 / ln_first_inward_mult;
  }
  const double nd = n;
  const double modal_k = nd * p_ddr.x[0];
  if (k + overflow_steps_lower_bound > modal_k) {
    // (ok, this duplicated code belongs in its own function...)
    const double min_mult_right = 1 + 0.03125 / (nmk * nmk);
    dd_real lastp_ddr = ddr_maked(1.0);
    dd_real right_sum_ddr = ddr_maked(0.0);
    while (1) {
      k += 1;
      lastp_ddr = ddr_mul(lastp_ddr, ddr_divd(ddr_muld(pdq_ddr, nmk), k));
      nmk -= 1;
      const double preaddp = right_sum_ddr.x[0];
      right_sum_ddr = ddr_add(right_sum_ddr, lastp_ddr);
      if (right_sum_ddr.x[0] <= preaddp * min_mult_right) {
        break;
      }
    }
    if (nmk > 0) {
      const double pdq = pdq_ddr.x[0];
      double lastp = lastp_ddr.x[0];
      double right_tail_sum = 0.0;
      while (1) {
        k += 1;
        lastp *= pdq * nmk / k;
        nmk -= 1;
        const double preaddp = right_tail_sum;
        right_tail_sum += lastp;
        if (right_tail_sum == preaddp) {
          break;
        }
      }
      right_sum_ddr = ddr_addd(right_sum_ddr, right_tail_sum);
    }
    k = obs_k;
    nmk = n - obs_k;
    lastp_ddr = ddr_maked(1.0);
    dd_real left_sum_ddr = lastp_ddr;
    if (k > 0) {
      const double min_mult_left = 1 + 0.0625 / (k * k);
      while (1) {
        nmk += 1;
        lastp_ddr = ddr_mul(lastp_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
        k -= 1;
        const double preaddp = left_sum_ddr.x[0];
        left_sum_ddr = ddr_add(left_sum_ddr, lastp_ddr);
        if (left_sum_ddr.x[0] <= preaddp * min_mult_left) {
          break;
        }
      }
      if (k > 0) {
        const double pdq = pdq_ddr.x[0];
        double lastp = lastp_ddr.x[0];
        double left_tail_sum = 0.0;
        while (1) {
          nmk += 1;
          lastp *= k / (pdq * nmk);
          k -= 1;
          const double preaddp = left_tail_sum;
          left_tail_sum += lastp;
          if (left_tail_sum == preaddp) {
            break;
          }
        }
        left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
      }
    }
    const dd_real prob_ddr = ddr_accurate_div(left_sum_ddr, ddr_add(left_sum_ddr, right_sum_ddr));
    if (!logp) {
      return prob_ddr.x[0];
    }
    return ddr_log(prob_ddr).x[0];
  }
  dd_real ln_prob_ddr = binom_ln_prob_internal(k, n, p_ddr);
  dd_real lastp_ddr = ddr_maked(1.0);
  dd_real left_sum_ddr = lastp_ddr;
  if (k > 0) {
    const double min_mult_left = 1 + 0.03125 / (k * k);
    while (1) {
      nmk += 1;
      lastp_ddr = ddr_mul(lastp_ddr, ddr_divd(ddr_muld(qdp_ddr, k), nmk));
      k -= 1;
      const double preaddp = left_sum_ddr.x[0];
      left_sum_ddr = ddr_add(left_sum_ddr, lastp_ddr);
      if (left_sum_ddr.x[0] <= preaddp * min_mult_left) {
        break;
      }
    }
    if (k > 0) {
      // Continue the calculation with ordinary precision.
      const double pdq = pdq_ddr.x[0];
      double lastp = lastp_ddr.x[0];
      double left_tail_sum = 0.0;
      while (1) {
        nmk += 1;
        lastp *= k / (pdq * nmk);
        k -= 1;
        const double preaddp = left_tail_sum;
        left_tail_sum += lastp;
        if (left_tail_sum == preaddp) {
          break;
        }
      }
      left_sum_ddr = ddr_addd(left_sum_ddr, left_tail_sum);
    }
  }
  ln_prob_ddr = ddr_add(ln_prob_ddr, ddr_log(left_sum_ddr));
  if (logp) {
    return ln_prob_ddr.x[0];
  }
  return ddr_exp(ln_prob_ddr).x[0];
}

// Assumes 0 <= n < 2^52, and should achieve quantile relative error < 2^{-54}
// unless n is well over 2^31.
// - The initial relative-likelihood should be accurate to at least ~60 bits
//   for n < 2^31.
// - We evaluate the inner part of the tailsum with full dd_real accuracy.  We
//   don't drop down to regular float64 accuracy until we know the additional
//   relative error introduced by doing so is < 2^{-56}.
// The goal is to support a higher-level qbinom() function where e.g.
//   qbinom(fractions.Fraction(1, 9), 2, fractions.Fraction(2, 3))
// can be trusted to be 0, and
//   qbinom(fractions.Fraction(2**53, (2**53 - 1) * 9), 2, fractions.Fraction(2, 3))
// can be trusted to be 1.
/*
int32_t Qbinom(dd_real ln_quantile_ddr, int64_t n, dd_real p_ddr) {
}
*/


#ifdef __cplusplus
}
#endif
