// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
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

#include "plink2_hwe.h"

#include <assert.h>
#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// - starting_lnprobv_ddr is expected to be initialized to
//     log(2^obs_hets / (obs_hets! obs_hom1! obs_hom2!))
//   or have x[0] initialized to DBL_MAX to indicate that the calculation
//   hasn't happened.  In the latter case, it may be set to the former value
//   if that is needed in the calculation.
//
// - *cmp_resultp is set to positive value if hets := obs_hets + 2*hom_decr has
//   higher likelihood than hets := obs_hets, 0 if identical likelihood, and
//   negative value if lower likelihood.
//
// - Error is returned iff malloc fails.
BoolErr HweCompare(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2, int32_t hom_decr, dd_real* starting_lnprobv_ddr_ptr, intptr_t* cmp_resultp, double* dbl_ptr) {
  // From e.g. the Wigginton paper, P(N_{AB}=n_{AB} | N, n_A) is
  //
  //      2^{n_{AB}} N! n_A! n_B!
  //   -----------------------------
  //   n_{AA}! n_{AB}! n_{BB}! (2N)!
  //
  // Thus, P(N_{AB}=obs_hets + 2*hom_decr) / P(N_{AB}=obs_hets) is
  //
  //       obs_hets!         obs_hom1! * obs_hom2!      2j
  //   ---------------- * -------------------------- * 2
  //   (obs_hets + 2j)!   (obs_hom1-j)!(obs_hom2-j)!
  //
  // where j=hom_decr.
  uint32_t numer_factorial_args[3];
  numer_factorial_args[0] = obs_hets;
  numer_factorial_args[1] = obs_hom1;
  numer_factorial_args[2] = obs_hom2;
  uint32_t denom_factorial_args[3];
  denom_factorial_args[0] = obs_hets + 2 * hom_decr;
  denom_factorial_args[1] = obs_hom1 - hom_decr;
  denom_factorial_args[2] = obs_hom2 - hom_decr;

  mp_limb_t* gmp_wkspace = nullptr;
  uintptr_t gmp_wkspace_limb_ct = 0;
  BoolErr reterr = CompareFactorialProducts(3, hom_decr * 2LL, obs_hets, numer_factorial_args, denom_factorial_args, starting_lnprobv_ddr_ptr, &gmp_wkspace, &gmp_wkspace_limb_ct, cmp_resultp, dbl_ptr);
  free_cond(gmp_wkspace);
  return reterr;
}

// obs_hets + obs_hom1 + obs_hom2 assumed to be <2^31.
BoolErr HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp, double* resultp) {
  // This function implements an exact SNP test of Hardy-Weinberg
  // Equilibrium as described in Wigginton, JE, Cutler, DJ, and
  // Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
  // Equilibrium. American Journal of Human Genetics. 76: 887 - 893.
  //
  // The original version was written by Jan Wigginton.
  //
  // This version was written by Christopher Chang.  It contains the following
  // improvements over the original SNPHWE():
  // 1. Proper handling of >64k genotypes.  Previously, there was a potential
  //    integer overflow.
  // 2. Detection and efficient handling of floating point overflow and
  //    underflow.  E.g. instead of summing a tail all the way down, the loop
  //    stops once the latest increment underflows the partial sum's 53-bit
  //    precision; this results in a large speedup when max heterozygote count
  //    >1k.
  // 3. No malloc() call in most cases: it's only necessary to keep track of a
  //    few partial sums.  (But see (6).)
  // 4. Support for the mid-p variant of this test.  See Graffelman J, Moreno V
  //    (2013) The mid p-value in exact tests for Hardy-Weinberg equilibrium.
  // 5. Log-p-value return (added Jan 2024).  A p-value of 1e-400 may be worth
  //    distinguishing from 1e-40000 in a biobank-scale dataset.
  // 6. Accurate handling of near-ties (added Mar 2026), using the QD and GMP
  //    libraries.  The QD library is used to efficiently perform computations
  //    on log-factorials with >60 bits of accuracy past the decimal point;
  //    this lets us jump from one tail to the other with negligible accuracy
  //    loss, and correctly resolve almost all near-ties without a full-blown
  //    rational expansion.  The GMP library supports the rational expansion
  //    when we need it; this can theoretically require large memory
  //    allocations, but that is very unlikely to come up in practice.
  //    (Note that we continue to allow the returned value to have a few bits
  //    of floating-point error.  But misclassification of a near-tie can
  //    result in a *large* relative error, so I've decided to go through the
  //    trouble of stamping that out despite its low analytical impact.)
  //
  // Note that the HweThreshLn() function is a lot more efficient for testing
  // against a p-value inclusion threshold.  HweLnP() should only be used if
  // you need the actual p-value.

  // Variables are mostly a mix of int32_ts and doubles, with a few
  // log-factorial-sum dd_reals (i.e. double-doubles).
  // int32 -> double casts are usually left implicit.
  // Naming conventions:
  // * obs_...: int32
  // * ..._ct: int32
  // * ..._ctd: double
  // * ..._ddr: dd_real
  // * Almost everything else is a double.
  int32_t obs_homc;
  int32_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  const int32_t rare_ct = 2 * obs_homr + obs_hets;
  if (rare_ct < 2) {
    *resultp = midp * (-kLn2);
    return 0;
  }
  // MAF: rare_ct / allele_ct
  // modal #hets:
  //   sample_ct * 2 * MAF * (1 - MAF)
  // = rare_ct * (1 - MAF)
  const int32_t sample_ct = obs_hom1 + obs_hom2 + obs_hets;
  // allele_ct can be >= 2^31.
  const double allele_ctd = S_CAST(double, sample_ct * 2LL);
  const double maf = rare_ct / allele_ctd;
  // possible todo: check whether this type of expression is worth rewriting as
  // e.g. prefer_fma(rare_ct, -maf, rare_ct)
  // 'c' in cmodal_nhet is for 'continuous', i.e. the mode if we extend the
  // likelihood function to all reals in [0, sample_ct].
  const double cmodal_nhet = rare_ct * (1 - maf);
  double hets = obs_hets;
  double homr = obs_homr;
  double homc = obs_homc;
  double lastp = 1;
  double tailp = 1 - midp * 0.5;
  if (hets > cmodal_nhet) {
    const double het_delta = hets - cmodal_nhet;
    // From Feb 2024 - Feb 2026, (except for the p=1 fast path) we always
    // computed starting log-likelihood.  Then, if it was high enough, we
    // proceeded with computing 1 - [sum of center likelihoods]; otherwise we
    // summed both tails in a manner that could handle likelihoods < DBL_MIN.
    //
    // Numerical stability of both branches was investigated in Mar 2026.
    // * The slightly-unstable 1 - [sum of center likelihoods] branch was
    //   reverted to the old cancellation-avoiding relative-likelihood
    //   algorithm, and the entrance condition changed to het_delta < 344.
    //   This doesn't risk center-likelihood overflow, since the
    //   relative-likelihood of the mode is loosely bounded above by
    //     ((172^172) / 172!)^4 ~= 5.3e+292
    //   which leaves enough headroom to accumulate the rest of the center-sum
    //   and multiply by e.g. allele_ct without overflowing.
    // * Log-factorial computations in the tail-jumping branch are now
    //   performed with "double-double" precision.  (Possible todo: benchmark
    //   float128 on x86_64.  But web search results imply QD is better.)
    if ((!midp) && (het_delta < 2.0)) {
      // Fast path for p=1.
      if (obs_hets * (obs_hets - 1LL) <= 4 * (obs_homc + 1LL) * (obs_homr + 1)) {
        *resultp = 0;
        return 0;
      }
    }
    // Iterate outward to floating-point precision limit.
    // No need for homr > 0 check, tailp == preaddp will trigger when we hit 0.
    // (hets, homr, and homc never accumulate any floating-point error since
    // they start as small integers and are only changed by adding/subtracting
    // 1 or 2.)
    while (1) {
      hets += 2;
      lastp *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
    }
    if (het_delta < 344.0) {
      // Jump back to starting contingency table, and iterate inward.
      lastp = 1;
      hets = obs_hets;
      homr = obs_homr;
      homc = obs_homc;
      dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
      double centerp = midp * 0.5;
      // No need for hets > 1 check, lastp checks do what we need.
      while (1) {
        homr += 1;
        homc += 1;
        lastp *= (hets * (hets - 1)) / (4 * homr * homc);
        hets -= 2;
        // Number of center tables is maximized with obs_hets - cmodal_nhet ~=
        // 344, obs_homr = 0, obs_homc and obs_hets both large.
        // Since 1 + 1/2 + ... + 1/172 < 1/173 + ... + 1/53000, we're limited
        // to ~53000 tables.  Each lastp update involves 4 operations which can
        // each introduce up to 0.5 ULP relative error under the default
        // rounding mode.
        if (lastp < 1 + 53000 * 2 * k2m52) {
          if (lastp <= 1 - 53000 * 2 * k2m52) {
            tailp += lastp;
            break;
          }
          // Near-tie.  True value of lastp can be greater than, equal to, or
          // less than 1.
          const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
          intptr_t cmp_result;
          if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
            return 1;
          }
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
      // No need for hets > 1 check, tailp == preaddp check does what we need
      // (even when hets is already -1 when entering the loop: in that case
      // lastp is 0).
      while (1) {
        homr += 1;
        homc += 1;
        lastp *= (hets * (hets - 1)) / (4 * homr * homc);
        hets -= 2;
        const double preaddp = tailp;
        tailp += lastp;
        if (tailp == preaddp) {
          break;
        }
      }
      *resultp = log(tailp / (tailp + centerp));
      return 0;
    }
    // starting_lnprobv_ddr is guaranteed to be negative for hets >= 4, and no
    // larger than ln(2) otherwise.
    const double c_minus_r = homc - homr;
    dd_real starting_lnprobv_ddr =
      ddr_sub(ddr_muld(_ddr_log2, obs_hets),
              ddr_add3_lfacts(obs_hets, obs_homr, obs_homc));
    // Now we want to jump near the other tail, without evaluating that many
    // contingency table log-likelihoods along the way.
    //
    // Each full log-likelihood evaluation requires 3 ddr_lfact() calls.  Since
    // they are now performed with extra precision, they require hundreds of
    // floating-point operations, so we want to limit ourselves to 1-2 full
    // evaluations most of the time.  (Possible todo: use lower-accuracy
    // Lfact() to jump around, followed by ddr_lfact() when exiting the loop.
    // Should be an easy performance win, but there's a complexity cost so I'll
    // wait until I see a scenario where this branch executes frequently...)
    //
    // The current heuristic starts by reflecting (obs_homr + homr) * 0.5
    // across the (continuous) mode, performing a full log-likelihood check at
    // an adjacent valid point.  Hopefully we find that we're in
    // (starting_lnprob - 62 * kLn2, starting_lnprob], so we're at or near a
    // table that actually contributes to the tail-sum.  (This window is chosen
    // to be wide enough to guarantee that at least one point falls inside when
    // sample_ct < 2^31.)
    //
    // If not, we jump again, using Newton's method.
    // If homr is too low (i.e. current log-likelihood is too high), when we
    // increase homr by 1, the likelihood gets multiplied by
    //   hets * (hets-1) / (4 * (homr+1) * (homc+1))
    // i.e. we're adding the logarithm of this value to the log-likelihood.
    // If homr is too high, when we decrease homr by 1, the likelihood gets
    // multiplied by
    //   4 * homr * homc / ((hets+2) * (hets+1))
    // We use the log of the first expression as the Newton's method f'(x) when
    // we're jumping to higher homr, and the negative-log of the second
    // expression when we're jumping to lower homr.
    // f''(x) is always negative, so we can aim for starting_lnprob instead of
    // the middle of the interval.

    // hets moves twice as fast as homr.  So if we add
    //   0.5 * (hets + obs_hets) - cmodal_nhet
    // to 0.5 * (homr + obs_homr), that reflects homr across the cmode.
    const double max_homr = S_CAST(double, rare_ct >> 1);
    {
      const double delta = 0.5 * (hets + obs_hets) - cmodal_nhet;
      homr = 0.5 * (homr + obs_homr) + delta;
      // Round up (to guarantee we've actually moved to the other side of the
      // cmode) and clamp.
      homr = ceil_smalleps_limit32(homr, max_homr);
    }
    // 'f' in lnprobf refers to fixed component of log-likelihood (unchanged
    // when rare_ct and sample_ct are held constant), 'v' refers to variable
    // component
    const dd_real lnprobf_ddr =
      ddr_sub(ddr_add3_lfacts(sample_ct, rare_ct, allele_ctd - rare_ct),
              ddr_lfact(allele_ctd));
    const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
    while (1) {
      hets = rare_ct - homr * 2;
      homc = homr + c_minus_r;
      const dd_real lnprobv_ddr =
        ddr_sub(ddr_muld(_ddr_log2, hets),
                ddr_add3_lfacts(hets, homr, homc));
      const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
      // Could tighten this threshold further; I haven't performed a careful
      // error analysis yet but CompareFactorialProducts() includes a plausible
      // assumption that 2^{-60} is safe.  But code is correct as long as we're
      // guaranteed to enter the "lastp < 2 - one_minus_scaled_eps" branch for
      // positive lnprob_diff.
      if (lnprob_diff >= k2m53) {
        if (homr == max_homr) {
          // All tables on this tail are larger than the starting table.  Exit.
          // (This is possible when obs_hom1 == obs_hom2 == 0.)
          *resultp = starting_lnprob + log(tailp);
          return 0;
        }
        const double ll_deriv = log(hets * (hets - 1) / (4 * (homr + 1) * (homc + 1)));
        // Round absolute value up, to guarantee that we make progress.
        // (lnprob_diff is positive and ll_deriv is negative.)
        // This may overshoot.  But the function is guaranteed to terminate
        // because we never overshoot (and we do always make progress on each
        // step) once we're on the other side.
        homr += ceil_smalleps(-lnprob_diff / ll_deriv);
        if (homr > max_homr) {
          homr = max_homr;
        }
      } else if (lnprob_diff > -62 * kLn2) {
        lastp = exp(lnprob_diff);
        break;
      } else {
        const double ll_deriv = log((hets + 2) * (hets + 1) / (4 * homr * homc));
        // Round down, to guarantee we don't overshoot.
        // We're guaranteed to make progress, since lnprob_diff <= -62 * log(2)
        // and sample_ct < 2^31.
        homr -= S_CAST(int64_t, lnprob_diff / ll_deriv);
      }
    }
    // Sum toward center, until lastp >= 1.  (No more risk of double-counting
    // the starting table, since we don't enter this branch at all unless the
    // starting table is >= 172 steps from the cmode.)
    //
    // lastp should be accurate to 3 ULP as we enter this loop (max 1.5 ULP
    // observed error from exp, tiny bit over 0.5 from lnprob_diff, we round up
    // all the way to 3 so we don't have to worry about "2 -
    // one_minus_scaled_eps" rounding behavior), so near-tie detection can use
    // a tight epsilon here.
    double one_minus_scaled_eps = 1 - 3 * k2m52;
    const double lastp_tail = lastp;
    const double homr_tail = homr;
    const double homc_tail = homc;
    const double hets_tail = hets;
    while (lastp <= one_minus_scaled_eps) {
      tailp += lastp;
      hets += 2;
      lastp *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lastp < 2 - one_minus_scaled_eps) {
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
      intptr_t cmp_result;
      if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
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
    homr = homr_tail;
    homc = homc_tail;
    hets = hets_tail;
    while (1) {
      homr += 1;
      homc += 1;
      lastp *= (hets * (hets - 1)) / (4 * homr * homc);
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
      hets -= 2;
    }
    *resultp = starting_lnprob + log(tailp);
    return 0;
  }
  // Same as above, just with directions flipped.
  const double het_delta = cmodal_nhet - hets;
  if ((!midp) && (het_delta < 2.0)) {
    // Fast path for p=1.
    if ((4LL * obs_homr) * obs_homc <= (obs_hets + 2LL) * (obs_hets + 1LL)) {
      *resultp = 0;
      return 0;
    }
  }
  // Iterate outward to floating-point precision limit.
  while (1) {
    homr += 1;
    homc += 1;
    lastp *= (hets * (hets - 1)) / (4 * homr * homc);
    hets -= 2;
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
  }
  if (het_delta < 344.0) {
    // Jump back to starting table, and iterate inward.
    lastp = 1;
    hets = obs_hets;
    homr = obs_homr;
    homc = obs_homc;
    dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
    double centerp = midp * 0.5;
    // No need for hets > 1 check, lastp checks do what we need.
    while (1) {
      hets += 2;
      lastp *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      // If we're 172 steps from the center, number of center tables is limited
      // to ~2*172 = 344, when obs_homr ~= obs_homc.
      if (lastp < 1 + 344 * 2 * k2m52) {
        if (lastp <= 1 - 344 * 2 * k2m52) {
          tailp += lastp;
          break;
        }
        // Near-tie.  True value of lastp can be greater than, equal to, or
        // less than 1.
        const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
        intptr_t cmp_result;
        if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
          return 1;
        }
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
      hets += 2;
      lastp *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
    }
    *resultp = log(tailp / (tailp + centerp));
    return 0;
  }
  const double c_minus_r = homc - homr;
  dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(_ddr_log2, obs_hets),
            ddr_add3_lfacts(obs_hets, obs_homr, obs_homc));
  // Jump to other tail.
  {
    const double delta = cmodal_nhet - 0.5 * (hets + obs_hets);
    homr = 0.5 * (homr + obs_homr) - delta;
    // Round down (to guarantee we've actually moved to the other side of the
    // cmode) and clamp.
    homr = S_CAST(int32_t, homr);
    if (homr < 0) {
      homr = 0;
    }
  }
#ifndef NDEBUG
  const double max_homr = S_CAST(double, rare_ct >> 1);
#endif
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_add3_lfacts(sample_ct, rare_ct, allele_ctd - rare_ct),
            ddr_lfact(allele_ctd));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  while (1) {
    hets = rare_ct - homr * 2;
    homc = homr + c_minus_r;
    const dd_real lnprobv_ddr =
      ddr_sub(ddr_muld(_ddr_log2, hets),
              ddr_add3_lfacts(hets, homr, homc));
    const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    if (lnprob_diff >= k2m53) {
      if (homr == 0) {
        // All tables on this tail have higher likelihood than the starting
        // table.  Exit.
        *resultp = starting_lnprob + log(tailp);
        return 0;
      }
      const double ll_deriv = log(4 * homr * homc / ((hets + 2) * (hets + 1)));
      homr -= ceil_smalleps(-lnprob_diff / ll_deriv);
      if (homr < 0) {
        homr = 0;
      }
    } else if (lnprob_diff > -62 * kLn2) {
      lastp = exp(lnprob_diff);
      break;
    } else {
      const double ll_deriv = log(4 * (homr + 1) * (homc + 1) / (hets * (hets - 1)));
      homr += S_CAST(int64_t, lnprob_diff / ll_deriv);
      assert(homr <= max_homr);
    }
  }
  // Sum toward center, until lastp >= 1.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double lastp_tail = lastp;
  const double homr_tail = homr;
  const double homc_tail = homc;
  const double hets_tail = hets;
  while (lastp <= one_minus_scaled_eps) {
    tailp += lastp;
    homr += 1;
    homc += 1;
    lastp *= hets * (hets - 1) / (4 * homr * homc);
    hets -= 2;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lastp < 2 - one_minus_scaled_eps) {
    const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
    intptr_t cmp_result;
    if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
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
  homr = homr_tail;
  homc = homc_tail;
  hets = hets_tail;
  while (1) {
    hets += 2;
    lastp *= 4 * homr * homc / (hets * (hets - 1));
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
    homr -= 1;
    homc -= 1;
  }
  *resultp = starting_lnprob + log(tailp);
  return 0;
}

// 2^{-83} bias to give plink 1.9-style exact tests maximum ability to
// determine tiny p-values.  (~2^{-53} is necessary to take advantage of
// denormalized small numbers, then allow tail sum to be up to 2^30.  ...okay,
// HweThresh[Midp]() is not responsible for denormal values of thresh, and
// plink2 now just flushes denormals to zero.  But configuring this
// constant to be compatible with them doesn't cost us anything.)
static const double kExactTestBias = k2m50 / (1LL << 33);

BoolErr HweThresh(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh, uint32_t* out_of_eqp) {
  // Threshold-test-only version of HweLnP() which is usually able to exit
  // from the calculation earlier.  Assumes DBL_MIN <= pval_thresh <= 1 (note
  // that some older versions of this function didn't handle pval_thresh=1
  // correctly, and plink2 still avoids calling this with pval_thresh=1).
  // Returns *out_of_eqp=0 if these counts are close enough to Hardy-Weinberg
  // equilibrium, 1 otherwise.
  //
  // BoolErr corresponds to malloc failure when resolving a near-tie.
  //
  // Suppose, for definiteness, that the number of observed hets is no less
  // than expectation.  (Same ideas apply for the other case.)  We proceed as
  // follows:
  // - Sum the *relative* likelihoods of more likely smaller het counts.
  // - Determine the minimum tail mass to pass the threshold.
  // - The majority of the time, the tail boundary elements are enough to pass
  //   the threshold; we never need to sum the remainder of the tails.
  // - And in the case of disequilibrium, we will often be able to immediately
  //   determine that the tail sum cannot possibly pass the threshold, just by
  //   looking at the tail boundary elements and using a geometric series to
  //   upper-bound the tail sums.
  // - Only when neither of these conditions hold do we start traveling down
  //   the tails.
  int32_t obs_homc;
  int32_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  const int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  if (!genotypes2) {
    *out_of_eqp = 0;
    return 0;
  }
  const int32_t rare_copies = 2 * obs_homr + obs_hets;
  double hets_t2 = obs_hets;  // tail 2
  double homr_t2 = obs_homr;
  double homc_t2 = obs_homc;

  double tailp1 = kExactTestBias;
  double centerp = 0;
  double lastp2 = tailp1;
  double tailp2 = 0;

  // const double center_div_tail_thresh = (1 - pval_thresh) / pval_thresh;
  // An initial upper bound on the tail sum is useful, since it lets us
  // report test failure before summing the entire center.  We use the
  // trivial bound of 1 + floor(rare_copies / 2): that's the total number
  // of possible het counts, and the relative probability for each count must
  // be <= 1 if it's in the tail.
  const double center_div_tail_thresh = (1 - pval_thresh) / pval_thresh;
  const double centerp_exit_thresh = u31tod(1 + (rare_copies >> 1)) * (center_div_tail_thresh * kExactTestBias);
  double scaled_one_plus_eps = kExactTestBias;

  // Expected het count:
  //   2 * rarefreq * (1 - rarefreq) * genotypes
  // = 2 * (rare_copies / (2 * genotypes)) * (1 - rarefreq) * genotypes
  // = rare_copies * (1 - (rare_copies / (2 * genotypes)))
  // = (rare_copies * (2 * genotypes - rare_copies)) / (2 * genotypes)
  //
  // The computational identity is
  //   P(nhets == n) := P(nhets == n+2) * (n+2) * (n+1) /
  //                    (4 * homr(n) * homc(n))
  // where homr() and homc() are the number of homozygous rares/commons needed
  // to maintain the same allele frequencies.
  // This probability is always decreasing when proceeding away from the
  // expected het count.
  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    // tail 1 = upper hets
    if (obs_hets < 2) {
      *out_of_eqp = 0;
      return 0;
    }

    // het_probs[hets] = 1
    // het_probs[hets - 2] = het_probs[hets] * hets * (hets - 1) / (4 * (homr + 1) * (homc + 1))
    do {
      homr_t2 += 1;
      homc_t2 += 1;
      lastp2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      hets_t2 -= 2;
      scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
      if (lastp2 < scaled_one_plus_eps) {
        if (lastp2 <= 2 * kExactTestBias - scaled_one_plus_eps) {
          tailp2 = lastp2;
          break;
        }
        double lastp = lastp2 * (1.0 / kExactTestBias);
        const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
        dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
        intptr_t cmp_result;
        if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        lastp2 = lastp * kExactTestBias;
        if (cmp_result <= 0) {
          tailp2 = lastp2;
          break;
        }
        // lastp can be of the form exp(lnprob_diff).  We are conservatively
        // assuming that exp() introduces up to 1.5 ULP error, and lnprob_diff
        // may be off by slightly over 0.5 ULP.
        scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
      }
      centerp += lastp2;
      if (centerp >= centerp_exit_thresh) {
        *out_of_eqp = 1;
        return 0;
      }
    } while (hets_t2 > 1);
    // hets_t2 guaranteed to be nonnegative on loop exit, so ratio
    // calculation works

    // This is NaN when pval_thresh=1, so we write the next if-condition to be
    // true on NaN.
    const double tailp_exit_thresh = centerp / center_div_tail_thresh;
    if (!(tailp1 + tailp2 < tailp_exit_thresh)) {
      *out_of_eqp = 0;
      return 0;
    }
    // c + cr + cr^2 + ... = c/(1-r), which is an upper bound for the tail sum
    const double ratio = (hets_t2 * (hets_t2 - 1)) / (4 * (homr_t2 + 1) * (homc_t2 + 1));
    const double tail2_ceil = tailp2 / (1 - ratio);
    double hets_t1 = obs_hets + 2;
    double homr_t1 = obs_homr;
    double homc_t1 = obs_homc;
    // ratio for the other tail
    double lastp1 = (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
    const double tail1_ceil = tailp1 / (1 - lastp1);
    if (tail1_ceil + tail2_ceil < tailp_exit_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
    lastp1 *= tailp1;
    tailp1 += lastp1;

    if (obs_homr > 1) {
      // het_probs[hets + 2] = het_probs[hets] * 4 * homr * homc / ((hets + 2) * (hets + 1))
      const double tailp1_exit_thresh = tailp_exit_thresh - tailp2;
      while (1) {
        hets_t1 += 2;
        homr_t1 -= 1;
        homc_t1 -= 1;
        lastp1 *= (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
        const double preaddp = tailp1;
        tailp1 += lastp1;
        if (tailp1 >= tailp1_exit_thresh) {
          *out_of_eqp = 0;
          return 0;
        }
        // homr_t1 == 1 check isn't necessary, but it provides a
        // noticeable speedup in my testing on real data.
        if ((tailp1 == preaddp) || (homr_t1 == 1)) {
          break;
        }
      }
    }
    if (tailp1 + tail2_ceil < tailp_exit_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
    const double tailp2_exit_thresh = tailp_exit_thresh - tailp1;
    while (1) {
      homr_t2 += 1;
      homc_t2 += 1;
      lastp2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      const double preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= tailp2_exit_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      if (tailp2 == preaddp) {
        *out_of_eqp = 1;
        return 0;
      }
      hets_t2 -= 2;
    }
  }
  // tail 1 = lower hets
  if (!obs_homr) {
    *out_of_eqp = 0;
    return 0;
  }
  do {
    hets_t2 += 2;
    lastp2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
    if (lastp2 < scaled_one_plus_eps) {
      if (lastp2 <= 2 * kExactTestBias - scaled_one_plus_eps) {
        tailp2 = lastp2;
        break;
      }
      double lastp = lastp2 * (1.0 / kExactTestBias);
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
      dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
      intptr_t cmp_result;
      if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
        return 1;
      }
      lastp2 = lastp * kExactTestBias;
      if (cmp_result <= 0) {
        tailp2 = lastp2;
        break;
      }
      scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
    }
    centerp += lastp2;
    if (centerp >= centerp_exit_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
  } while (homr_t2 > 0);
  // homr_t2 guaranteed to be nonnegative on loop exit, so ratio
  // calculation works

  const double tailp_exit_thresh = centerp / center_div_tail_thresh;
  if (!(tailp1 + tailp2 < tailp_exit_thresh)) {
    *out_of_eqp = 0;
    return 0;
  }
  const double ratio = (4 * homr_t2 * homc_t2) / ((hets_t2 + 2) * (hets_t2 + 1));
  const double tail2_ceil = tailp2 / (1 - ratio);
  double hets_t1 = obs_hets;
  double homr_t1 = obs_homr + 1;
  double homc_t1 = obs_homc + 1;
  double lastp1 = (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
  const double tail1_ceil = tailp1 / (1 - lastp1);
  if (tail1_ceil + tail2_ceil < tailp_exit_thresh) {
    *out_of_eqp = 1;
    return 0;
  }
  lastp1 *= tailp1;
  tailp1 += lastp1;

  if (obs_hets >= 4) {
    const double tailp1_exit_thresh = tailp_exit_thresh - tailp2;
    while (1) {
      hets_t1 -= 2;
      homr_t1 += 1;
      homc_t1 += 1;
      lastp1 *= (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
      const double preaddp = tailp1;
      tailp1 += lastp1;
      if (tailp1 >= tailp1_exit_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      if ((tailp1 == preaddp) || (hets_t1 < 4)) {
        break;
      }
    }
  }
  if (tailp1 + tail2_ceil < tailp_exit_thresh) {
    *out_of_eqp = 1;
    return 0;
  }
  const double tailp2_exit_thresh = tailp_exit_thresh - tailp1;
  while (1) {
    hets_t2 += 2;
    lastp2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    const double preaddp = tailp2;
    tailp2 += lastp2;
    if (tailp2 >= tailp2_exit_thresh) {
      *out_of_eqp = 0;
      return 0;
    }
    if (tailp2 == preaddp) {
      *out_of_eqp = 1;
      return 0;
    }
  }
}

BoolErr HweThreshMidp(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh, uint32_t* out_of_eqp) {
  // Mid-p version of HweThresh().  (There are enough fiddly differences that I
  // think it's better for this to be a separate function.)  Assumes
  // DBL_MIN <= pval_thresh < 0.5.
  int32_t obs_homc;
  int32_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  const int64_t genotypes2 = (obs_hets + obs_homc + obs_homr) * 2LL;
  if (!genotypes2) {
    *out_of_eqp = 0;
    return 0;
  }
  int32_t rare_copies = 2 * obs_homr + obs_hets;
  double hets_t2 = obs_hets;  // tail 2
  double homr_t2 = obs_homr;
  double homc_t2 = obs_homc;
  double tailp1 = kExactTestBias * 0.5;
  double centerp = tailp1;
  double lastp2 = kExactTestBias;
  double tailp2 = 0;
  const double center_div_tail_thresh = (1 - pval_thresh) / pval_thresh;
  const double centerp_exit_thresh = u31tod(1 + (rare_copies >> 1)) * (center_div_tail_thresh * kExactTestBias);
  double scaled_one_plus_eps = kExactTestBias;
  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    if (obs_hets < 2) {
      *out_of_eqp = 0;
      return 0;
    }
    do {
      homr_t2 += 1;
      homc_t2 += 1;
      lastp2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      hets_t2 -= 2;
      scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
      if (lastp2 < scaled_one_plus_eps) {
        if (lastp2 <= 2 * kExactTestBias - scaled_one_plus_eps) {
          tailp2 = lastp2;
          break;
        }
        double lastp = lastp2 * (1.0 / kExactTestBias);
        const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
        dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
        intptr_t cmp_result;
        if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        lastp2 = lastp * kExactTestBias;
        if (cmp_result <= 0) {
          if (cmp_result == 0) {
            tailp2 = tailp1;
            centerp += tailp1;
          } else {
            tailp2 = lastp2;
          }
          break;
        }
        scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
      }
      centerp += lastp2;
      if (centerp >= centerp_exit_thresh) {
        *out_of_eqp = 1;
        return 0;
      }
    } while (hets_t2 > 1);

    const double tailp_exit_thresh = centerp / center_div_tail_thresh;
    if (!(tailp1 + tailp2 < tailp_exit_thresh)) {
      *out_of_eqp = 0;
      return 0;
    }
    const double ratio = (hets_t2 * (hets_t2 - 1)) / (4 * (homr_t2 + 1) * (homc_t2 + 1));
    // this needs to work in both the tie and no-tie cases
    const double tail2_ceil = prefer_fma(lastp2, ratio / (1 - ratio), tailp2);
    double hets_t1 = obs_hets + 2;
    double homr_t1 = obs_homr;
    double homc_t1 = obs_homc;
    double lastp1 = (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
    // always a tie here
    const double tail1_ceil = (tailp1 * 2) / (1 - lastp1) - tailp1;
    if (tail1_ceil + tail2_ceil < tailp_exit_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
    lastp1 *= tailp1 * 2;
    tailp1 += lastp1;

    if (obs_homr > 1) {
      const double tailp1_exit_thresh = tailp_exit_thresh - tailp2;
      while (1) {
        hets_t1 += 2;
        homr_t1 -= 1;
        homc_t1 -= 1;
        lastp1 *= (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
        const double preaddp = tailp1;
        tailp1 += lastp1;
        if (tailp1 >= tailp1_exit_thresh) {
          *out_of_eqp = 0;
          return 0;
        }
        if ((tailp1 == preaddp) || (homr_t1 == 1)) {
          break;
        }
      }
    }
    if (tailp1 + tail2_ceil < tailp_exit_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
    const double tailp2_exit_thresh = tailp_exit_thresh - tailp1;
    while (1) {
      homr_t2 += 1;
      homc_t2 += 1;
      lastp2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      const double preaddp = tailp2;
      tailp2 += lastp2;
      if (tailp2 >= tailp2_exit_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      if (tailp2 == preaddp) {
        *out_of_eqp = 1;
        return 0;
      }
      hets_t2 -= 2;
    }
  }
  if (!obs_homr) {
    *out_of_eqp = 0;
    return 0;
  }
  do {
    hets_t2 += 2;
    lastp2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
    if (lastp2 < scaled_one_plus_eps) {
      if (lastp2 <= 2 - scaled_one_plus_eps) {
        tailp2 = lastp2;
        break;
      }
      double lastp = lastp2 * (1.0 / kExactTestBias);
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
      dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
      intptr_t cmp_result;
      if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
        return 1;
      }
      lastp2 = lastp * kExactTestBias;
      if (cmp_result <= 0) {
        if (cmp_result == 0) {
          tailp2 = tailp1;
          centerp += tailp1;
        } else {
          tailp2 = lastp2;
        }
        break;
      }
      scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
    }
    centerp += lastp2;
    if (centerp >= centerp_exit_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
  } while (homr_t2 > 0);

  const double tailp_exit_thresh = centerp / center_div_tail_thresh;
  if (!(tailp1 + tailp2 < tailp_exit_thresh)) {
    *out_of_eqp = 0;
    return 0;
  }
  const double ratio = (4 * homr_t2 * homc_t2) / ((hets_t2 + 2) * (hets_t2 + 1));
  const double tail2_ceil = prefer_fma(lastp2, ratio / (1 - ratio), tailp2);
  double hets_t1 = obs_hets;
  double homr_t1 = obs_homr + 1;
  double homc_t1 = obs_homc + 1;
  double lastp1 = (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
  const double tail1_ceil = (2 * tailp1) / (1 - lastp1) - tailp1;
  lastp1 *= 2 * tailp1;
  tailp1 += lastp1;

  if (tail1_ceil + tail2_ceil < tailp_exit_thresh) {
    *out_of_eqp = 1;
    return 0;
  }
  if (obs_hets >= 4) {
    const double tailp1_exit_thresh = tailp_exit_thresh - tailp2;
    while (1) {
      hets_t1 -= 2;
      homr_t1 += 1;
      homc_t1 += 1;
      lastp1 *= (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
      const double preaddp = tailp1;
      tailp1 += lastp1;
      if (tailp1 >= tailp1_exit_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      if ((tailp1 == preaddp) || (hets_t1 < 4)) {
        break;
      }
    }
  }
  if (tailp1 + tail2_ceil < tailp_exit_thresh) {
    *out_of_eqp = 1;
    return 0;
  }
  const double tailp2_exit_thresh = tailp_exit_thresh - tailp1;
  while (1) {
    hets_t2 += 2;
    lastp2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    const double preaddp = tailp2;
    tailp2 += lastp2;
    if (tailp2 >= tailp2_exit_thresh) {
      *out_of_eqp = 0;
      return 0;
    }
    if (tailp2 == preaddp) {
      *out_of_eqp = 1;
      return 0;
    }
  }
}

BoolErr HweThreshLnMain(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp, double ln_thresh, uint32_t* out_of_eqp) {
  assert(ln_thresh < -708);
  // Threshold-test-only version of HweLnP() which is usually able to exit
  // from the calculation earlier.  Returns *out_of_eqp=0 if these counts are
  // close enough to Hardy-Weinberg equilibrium, 1 otherwise.
  //
  // Assumes ln_thresh < -708, otherwise the |curr_hets - cmodal_nhet| < 344
  // early-exit doesn't work.  (possible todo: make the 344 constant a function
  // of ln_thresh, and set it closer to optimally.)
  //
  // Caller is responsible for including a tolerance in ln_thresh when
  // appropriate.
  int32_t obs_homc;
  int32_t obs_homr;
  if (obs_hom1 < obs_hom2) {
    obs_homc = obs_hom2;
    obs_homr = obs_hom1;
  } else {
    obs_homc = obs_hom1;
    obs_homr = obs_hom2;
  }
  const int32_t rare_ct = 2 * obs_homr + obs_hets;
  const int32_t sample_ct = obs_hom1 + obs_hom2 + obs_hets;
  const double allele_ctd = S_CAST(double, sample_ct * 2LL);
  const double maf = rare_ct / allele_ctd;
  const double cmodal_nhet = rare_ct * (1 - maf);
  double hets = obs_hets;
  // Change this to "rare_ct < 2" if ln_thresh restriction is being loosened
  // (to e.g. compare results against HweThresh()).
  if (fabs(hets - cmodal_nhet) < 344.0) {
    *out_of_eqp = 0;
    return 0;
  }

  // 1. Compute log-likelihood of starting table.  This may be high enough on
  //    its own to immediately return 0.
  //    If likelihood is lower than threshold / <total # of tables>, we can
  //    immediately return 1.
  // 2. Determine tailsum we must hit to return 0.
  // 3. The rest follows HweLnP(), except with an extra geometric-series-based
  //    early-exit attempt near the end.
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_add3_lfacts(sample_ct, rare_ct, allele_ctd - rare_ct),
            ddr_lfact(allele_ctd));
  double homr = obs_homr;
  double homc = obs_homc;
  dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(_ddr_log2, hets),
            ddr_add3_lfacts(hets, homr, homc));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  if (ln_thresh <= starting_lnprob - midp * kLn2) {
    *out_of_eqp = 0;
    return 0;
  }
  const double max_homr = S_CAST(double, rare_ct >> 1);
  if (ln_thresh > starting_lnprob + log(max_homr + 1 - midp * 0.5)) {
    *out_of_eqp = 1;
    return 0;
  }

  const double c_minus_r = homc - homr;
  // This should be in (0.5, 2^30].
  const double tail_thresh = exp(ln_thresh - starting_lnprob);
  double tailp = 1 - midp * 0.5;
  double lastp = 1;
  if (hets > cmodal_nhet) {
    // No center-sum (or p=1) code path, since it doesn't make sense to choose
    // a HWE threshold that makes these relevant; we should have already
    // returned 0.
    // (So we can't assume obs_hets >= 2.)
    while (1) {
      hets += 2;
      lastp *= (4 * homr * homc) / (hets * (hets - 1));
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
      if (tailp >= tail_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      homr -= 1;
      homc -= 1;
    }
    {
      const double delta = 0.5 * (hets + obs_hets) - cmodal_nhet;
      homr = 0.5 * (homr + obs_homr) + delta;
      homr = ceil_smalleps_limit32(homr, max_homr);
    }
    while (1) {
      hets = rare_ct - homr * 2;
      homc = homr + c_minus_r;
      const dd_real lnprobv_ddr =
        ddr_sub(ddr_muld(_ddr_log2, hets),
                ddr_add3_lfacts(hets, homr, homc));
      const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
      if (lnprob_diff >= k2m53) {
        if (homr == max_homr) {
          *out_of_eqp = 1;
          return 0;
        }
        const double ll_deriv = log(hets * (hets - 1) / (4 * (homr + 1) * (homc + 1)));
        homr += ceil_smalleps(-lnprob_diff / ll_deriv);
        if (homr > max_homr) {
          homr = max_homr;
        }
      } else if (lnprob_diff > -62 * kLn2) {
        lastp = exp(lnprob_diff);
        break;
      } else {
        const double ll_deriv = log((hets + 2) * (hets + 1) / (4 * homr * homc));
        homr -= S_CAST(int64_t, lnprob_diff / ll_deriv);
        assert(homr >= 0);
      }
    }
    double one_minus_scaled_eps = 1 - 3 * k2m52;
    const double lastp_tail = lastp;
    const double homr_tail = homr;
    const double homc_tail = homc;
    const double hets_tail = hets;
    while (lastp <= one_minus_scaled_eps) {
      tailp += lastp;
      if (tailp >= tail_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      hets += 2;
      lastp *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lastp < 2 - one_minus_scaled_eps) {
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
      intptr_t cmp_result;
      if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
        return 1;
      }
      if (cmp_result <= 0) {
        if (cmp_result == 0) {
          tailp += 1 - midp * 0.5;
        } else {
          tailp += lastp;
        }
        if (tailp >= tail_thresh) {
          *out_of_eqp = 0;
          return 0;
        }
      }
    }
    // We're down to one tail that can be tightly bounded by a geometric
    // series.  (ratio is always decreasing)
    // c + cr + cr^2 + ... = c/(1-r)
    lastp = lastp_tail;
    homr = homr_tail;
    homc = homc_tail;
    hets = hets_tail;

    homr += 1;
    homc += 1;
    const double cur_ratio = (hets * (hets - 1)) / (4 * homr * homc);
    lastp *= cur_ratio;
    const double remaining_ceil = lastp_tail / (1 - cur_ratio);
    if (tailp + remaining_ceil < tail_thresh) {
      *out_of_eqp = 1;
      return 0;
    }
    while (1) {
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        *out_of_eqp = 1;
        return 0;
      }
      if (tailp >= tail_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
      hets -= 2;
      homr += 1;
      homc += 1;
      lastp *= (hets * (hets - 1)) / (4 * homr * homc);
    }
    // unreachable
  }
  while (1) {
    homr += 1;
    homc += 1;
    lastp *= hets * (hets - 1) / (4 * homr * homc);
    hets -= 2;
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
    if (tailp >= tail_thresh) {
      *out_of_eqp = 0;
      return 0;
    }
  }
  {
    const double delta = cmodal_nhet - 0.5 * (hets + obs_hets);
    homr = 0.5 * (homr + obs_homr) - delta;
    homr = S_CAST(int32_t, homr);
    if (homr < 0) {
      homr = 0;
    }
  }
  while (1) {
    hets = rare_ct - homr * 2;
    homc = homr + c_minus_r;
    const dd_real lnprobv_ddr =
      ddr_sub(ddr_muld(_ddr_log2, hets),
              ddr_add3_lfacts(hets, homr, homc));
    const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    if (lnprob_diff >= k2m53) {
      if (homr <= 0) {
        *out_of_eqp = 1;
        return 0;
      }
      const double ll_deriv = log(4 * homr * homc / ((hets + 2) * (hets + 1)));
      homr -= ceil_smalleps(-lnprob_diff / ll_deriv);
      if (homr < 0) {
        homr = 0;
      }
    } else if (lnprob_diff > -62 * kLn2) {
      lastp = exp(lnprob_diff);
      break;
    } else {
      const double ll_deriv = log(4 * (homr + 1) * (homc + 1) / (hets * (hets - 1)));
      homr += S_CAST(int64_t, lnprob_diff / ll_deriv);
      assert(homr <= max_homr);
    }
  }
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double lastp_tail = lastp;
  const double homr_tail = homr;
  const double homc_tail = homc;
  const double hets_tail = hets;
  while (lastp <= one_minus_scaled_eps) {
    tailp += lastp;
    if (tailp >= tail_thresh) {
      *out_of_eqp = 0;
      return 0;
    }
    homr += 1;
    homc += 1;
    lastp *= hets * (hets - 1) / (4 * homr * homc);
    hets -= 2;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lastp < 2 - one_minus_scaled_eps) {
    const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
    intptr_t cmp_result;
    if (unlikely(HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
      return 1;
    }
    if (cmp_result <= 0) {
      if (cmp_result == 0) {
        tailp += 1 - midp * 0.5;
      } else {
        tailp += lastp;
      }
      if (tailp >= tail_thresh) {
        *out_of_eqp = 0;
        return 0;
      }
    }
  }
  lastp = lastp_tail;
  homr = homr_tail;
  homc = homc_tail;
  hets = hets_tail;

  hets += 2;
  const double cur_ratio = 4 * homr * homc / (hets * (hets - 1));
  lastp *= cur_ratio;
  const double remaining_ceil = lastp / (1 - cur_ratio);
  if (tailp + remaining_ceil < tail_thresh) {
    *out_of_eqp = 1;
    return 0;
  }
  while (1) {
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      *out_of_eqp = 1;
      return 0;
    }
    if (tailp >= tail_thresh) {
      *out_of_eqp = 0;
      return 0;
    }
    homr -= 1;
    homc -= 1;
    hets += 2;
    lastp *= 4 * homr * homc / (hets * (hets - 1));
  }
}

#ifdef __cplusplus
}
#endif
