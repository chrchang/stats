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

// - If non-null, *starting_lnprobv_ddr_ptr is expected to be initialized to
//     log(2^obs_hets / (obs_hets! obs_hom1! obs_hom2!)).
//
// - Returns positive value if hets := obs_hets + 2*hom_decr has higher
//   probability than hets := obs_hets, 0 if identical probability, and
//   negative value if lower probability.
intptr_t HweCompare(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2, int32_t hom_decr, dd_real* starting_lnprobv_ddr_ptr, double* dbl_ptr) {
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
  uint64_t numer_factorial_args[3];
  numer_factorial_args[0] = obs_hets;
  numer_factorial_args[1] = obs_hom1;
  numer_factorial_args[2] = obs_hom2;
  uint64_t denom_factorial_args[3];
  denom_factorial_args[0] = obs_hets + 2 * hom_decr;
  denom_factorial_args[1] = obs_hom1 - hom_decr;
  denom_factorial_args[2] = obs_hom2 - hom_decr;
  td_real ln_odds_ratio_tdr = _tdr_log2;
  td_real starting_lnprobv_tdr;
  if (starting_lnprobv_ddr_ptr == nullptr) {
    starting_lnprobv_tdr = tdr_make1(DBL_MAX);
  } else {
    starting_lnprobv_tdr = tdr_make(starting_lnprobv_ddr_ptr->x[0], starting_lnprobv_ddr_ptr->x[1], DBL_MAX);
  }
  return CompareFactorialProducts(3, tdr_make1(2.0), hom_decr * 2LL, obs_hets, numer_factorial_args, denom_factorial_args, &starting_lnprobv_tdr, &ln_odds_ratio_tdr, dbl_ptr);
}

// obs_hets + obs_hom1 + obs_hom2 assumed to be <2^31.
double HweLnP(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp) {
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
  // 6. Highly accurate handling of near-ties (added Mar-May 2026), using the
  //    QD library.  The QD library is used to efficiently perform computations
  //    on log-factorials with >60 bits of accuracy past the decimal point;
  //    this lets us jump from one tail to the other with negligible accuracy
  //    loss, and correctly resolve almost all near-ties.
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
    return midp * (-kLn2);
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
  double lik = 1;
  double tail_sum = 1 - midp * 0.5;
  if (hets > cmodal_nhet) {
    const double het_delta = hets - cmodal_nhet;
    // From Feb 2024 - Feb 2026, (except for the p=1 fast path) we always
    // computed starting log-probability.  Then, if it was high enough, we
    // proceeded with computing 1 - [sum of center probabilities]; otherwise we
    // summed both tails in a manner that could handle probabilities < DBL_MIN.
    //
    // Numerical stability of both branches was investigated in Mar 2026.
    // * The slightly-unstable 1 - [sum of center probabilities] branch was
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
        return 0.0;
      }
    }
    // Iterate outward to floating-point precision limit.
    // No need for homr > 0 check, tail_sum == preadd will trigger when we hit
    // 0.
    // (hets, homr, and homc never accumulate any floating-point error since
    // they start as small integers and are only changed by adding/subtracting
    // 1 or 2.)
    while (1) {
      hets += 2;
      lik *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
    }
    if (het_delta < 344.0) {
      // Jump back to starting contingency table, and iterate inward.
      lik = 1;
      hets = obs_hets;
      homr = obs_homr;
      homc = obs_homc;
      double center_sum = midp * 0.5;
      // No need for hets > 1 check, lik checks do what we need.
      while (1) {
        homr += 1;
        homc += 1;
        lik *= (hets * (hets - 1)) / (4 * homr * homc);
        hets -= 2;
        // Number of center tables is maximized with obs_hets - cmodal_nhet ~=
        // 344, obs_homr = 0, obs_homc and obs_hets both large.
        // Since 1 + 1/2 + ... + 1/172 < 1/173 + ... + 1/53000, we're limited
        // to ~53000 tables.  Each lik update involves 4 operations which can
        // each introduce up to 0.5 ULP relative error under the default
        // rounding mode.  One ULP <= 2^{-52};
        //   (1 + 2^{-52})^k < 1 + (k+1)2^{-52}
        // for k < tens of millions, so in most cases we can safely add
        // error-bounds together as long as we add an extra 2^{-52} at the
        // beginning.
        if (lik < 1 + 53000 * 2 * k2m52) {
          if (lik <= 1 - 53000 * 2 * k2m52) {
            tail_sum += lik;
            break;
          }
          // Near-tie.  True value of lik can be greater than, equal to, or
          // less than 1.
          const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
          const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, nullptr, &lik);
          if (cmp_result <= 0) {
            tail_sum += lik;
            if (midp && (cmp_result == 0)) {
              tail_sum -= 0.5;
              center_sum += 0.5;
            }
            break;
          }
        }
        center_sum += lik;
      }
      // Continue down tail to floating-point precision limit.
      // No need for hets > 1 check, tail_sum == preadd check does what we need
      // (even when hets is already -1 when entering the loop: in that case lik
      // is 0).
      while (1) {
        homr += 1;
        homc += 1;
        lik *= (hets * (hets - 1)) / (4 * homr * homc);
        hets -= 2;
        const double preadd = tail_sum;
        tail_sum += lik;
        if (tail_sum == preadd) {
          break;
        }
      }
      return log(tail_sum / (tail_sum + center_sum));
    }
    // starting_lnprobv_ddr is guaranteed to be negative for hets >= 4, and no
    // larger than ln(2) otherwise.
    const double c_minus_r = homc - homr;
    dd_real starting_lnprobv_ddr =
      ddr_sub(ddr_muld(_ddr_log2, obs_hets),
              ddr_add3_lfacts(obs_homr, obs_hets, obs_homc));
    // Now we want to jump near the other tail, without evaluating that many
    // contingency table log-probabilities along the way.
    //
    // Each full log-probability evaluation requires 3 ddr_lfact() calls.
    // Since they are now performed with extra precision, they require hundreds
    // of floating-point operations, so we want to limit ourselves to 1-2 full
    // evaluations most of the time.  (Possible todo: use lower-accuracy
    // Lfact() to jump around, followed by ddr_lfact() when exiting the loop.
    // Should be an easy performance win, but there's a complexity cost so I'll
    // wait until I see a scenario where this branch executes frequently...)
    //
    // The current heuristic starts by reflecting (obs_homr + homr) * 0.5
    // across the (continuous) mode, performing a full log-probability check at
    // an adjacent valid point.  Hopefully we find that we're in
    // (starting_lnprob - 62 * kLn2, starting_lnprob], so we're at or near a
    // table that actually contributes to the tail-sum float64.  (This window
    // is chosen to be wide enough to guarantee that at least one point falls
    // inside when sample_ct < 2^31.)
    //
    // If not, we jump again, using Newton's method.
    // If homr is too low (i.e. current log-probability is too high), when we
    // increase homr by 1, the probability gets multiplied by
    //   hets * (hets-1) / (4 * (homr+1) * (homc+1))
    // i.e. we're adding the logarithm of this value to the log-probability.
    // If homr is too high, when we decrease homr by 1, the probability gets
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
      homr = ceil_limit(homr, max_homr);
    }
    // 'f' in lnprobf refers to fixed component of log-probability (unchanged
    // when rare_ct and sample_ct are held constant), 'v' refers to variable
    // component
    const dd_real lnprobf_ddr =
      ddr_sub(ddr_add3_lfacts(rare_ct, sample_ct, allele_ctd - rare_ct),
              ddr_lfact(allele_ctd));
    const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
    while (1) {
      hets = rare_ct - homr * 2;
      homc = homr + c_minus_r;
      const dd_real lnprobv_ddr =
        ddr_sub(ddr_muld(_ddr_log2, hets),
                ddr_add3_lfacts(homr, hets, homc));
      const double lnprobv_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
      // Could tighten this threshold further.  But code is correct as long as
      // we're guaranteed to enter the "lik < 2 - one_minus_scaled_eps" branch
      // for positive lnprobv_diff.
      if (lnprobv_diff >= k2m60) {
        if (homr == max_homr) {
          // All tables on this tail are larger than the starting table.  Exit.
          // (This is possible when obs_hom1 == obs_hom2 == 0.)
          return starting_lnprob + log(tail_sum);
        }
        const double lnprobv_deriv = log(hets * (hets - 1) / (4 * (homr + 1) * (homc + 1)));
        // Round absolute value up, to guarantee that we make progress.
        // (lnprobv_diff is positive and lnprobv_deriv is negative.)
        // This may overshoot.  But the function is guaranteed to terminate
        // because we never overshoot (and we do always make progress on each
        // step) once we're on the other side.
        homr += ceil(-lnprobv_diff / lnprobv_deriv);
        if (homr > max_homr) {
          homr = max_homr;
        }
      } else if (lnprobv_diff > -62 * kLn2) {
        lik = exp(lnprobv_diff);
        break;
      } else {
        const double lnprobv_deriv = log((hets + 2) * (hets + 1) / (4 * homr * homc));
        // Round down, to guarantee we don't overshoot.
        // We're guaranteed to make progress, since lnprobv_diff <=
        // -62 * log(2) and sample_ct < 2^31.
        homr -= S_CAST(int64_t, lnprobv_diff / lnprobv_deriv);
      }
    }
    // Sum toward center, until lik >= 1.  (No more risk of double-counting the
    // starting table, since we don't enter this branch at all unless the
    // starting table is >= 172 steps from the cmode.)
    //
    // lik should be accurate to 3 ULP as we enter this loop (max 1.5 ULP
    // observed error from exp, tiny bit over 0.5 from lnprobv_diff, we round
    // up all the way to 3 so we don't have to worry about "2 -
    // one_minus_scaled_eps" rounding behavior), so near-tie detection can use
    // a tight epsilon here.
    double one_minus_scaled_eps = 1 - 3 * k2m52;
    // Save where we're starting on this tail, which isn't necessarily on the
    // boundary.  We sum inward until relative-likelihood > 1, then we jump
    // back to tailenter_{hets,homr,homc} and sum outward.
    const double tailenter_lik = lik;
    const double tailenter_homr = homr;
    const double tailenter_homc = homc;
    const double tailenter_hets = hets;
    while (lik <= one_minus_scaled_eps) {
      tail_sum += lik;
      hets += 2;
      lik *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lik < 2 - one_minus_scaled_eps) {
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
      const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &lik);
      if (cmp_result <= 0) {
        tail_sum += lik;
        if (midp && (cmp_result == 0)) {
          tail_sum -= 0.5;
        }
      }
    }
    // Sum away from center, until sums stop changing.
    lik = tailenter_lik;
    homr = tailenter_homr;
    homc = tailenter_homc;
    hets = tailenter_hets;
    while (1) {
      homr += 1;
      homc += 1;
      lik *= (hets * (hets - 1)) / (4 * homr * homc);
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
      hets -= 2;
    }
    return starting_lnprob + log(tail_sum);
  }
  // Same as above, just with directions flipped.
  const double het_delta = cmodal_nhet - hets;
  if ((!midp) && (het_delta < 2.0)) {
    // Fast path for p=1.
    if ((4LL * obs_homr) * obs_homc <= (obs_hets + 2LL) * (obs_hets + 1LL)) {
      return 0.0;
    }
  }
  // Iterate outward to floating-point precision limit.
  while (1) {
    homr += 1;
    homc += 1;
    lik *= (hets * (hets - 1)) / (4 * homr * homc);
    hets -= 2;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
  }
  if (het_delta < 344.0) {
    // Jump back to starting table, and iterate inward.
    lik = 1;
    hets = obs_hets;
    homr = obs_homr;
    homc = obs_homc;
    double center_sum = midp * 0.5;
    // No need for hets > 1 check, lik checks do what we need.
    while (1) {
      hets += 2;
      lik *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      // If we're 172 steps from the center, number of center tables is limited
      // to ~2*172 = 344, when obs_homr ~= obs_homc.
      if (lik < 1 + (344 * 2 + 1) * k2m52) {
        if (lik <= 1 - (344 * 2 + 1) * k2m52) {
          tail_sum += lik;
          break;
        }
        // Near-tie.  True value of lik can be greater than, equal to, or less
        // than 1.
        const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
        const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, nullptr, &lik);
        if (cmp_result <= 0) {
          tail_sum += lik;
          if (midp && (cmp_result == 0)) {
            tail_sum -= 0.5;
            center_sum += 0.5;
          }
          break;
        }
      }
      center_sum += lik;
    }
    // Continue down tail to floating-point precision limit.
    while (1) {
      hets += 2;
      lik *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
    }
    return log(tail_sum / (tail_sum + center_sum));
  }
  const double c_minus_r = homc - homr;
  dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(_ddr_log2, obs_hets),
            ddr_add3_lfacts(obs_homr, obs_hets, obs_homc));
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
    ddr_sub(ddr_add3_lfacts(rare_ct, sample_ct, allele_ctd - rare_ct),
            ddr_lfact(allele_ctd));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  while (1) {
    hets = rare_ct - homr * 2;
    homc = homr + c_minus_r;
    const dd_real lnprobv_ddr =
      ddr_sub(ddr_muld(_ddr_log2, hets),
              ddr_add3_lfacts(homr, hets, homc));
    const double lnprobv_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    if (lnprobv_diff >= k2m60) {
      if (homr == 0) {
        // All tables on this tail have higher probability than the starting
        // table.  Exit.
        return starting_lnprob + log(tail_sum);
      }
      const double lnprobv_deriv = log(4 * homr * homc / ((hets + 2) * (hets + 1)));
      homr -= ceil(-lnprobv_diff / lnprobv_deriv);
      if (homr < 0) {
        homr = 0;
      }
    } else if (lnprobv_diff > -62 * kLn2) {
      lik = exp(lnprobv_diff);
      break;
    } else {
      const double lnprobv_deriv = log(4 * (homr + 1) * (homc + 1) / (hets * (hets - 1)));
      homr += S_CAST(int64_t, lnprobv_diff / lnprobv_deriv);
      assert(homr <= max_homr);
    }
  }
  // Sum toward center, until lik >= 1.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double tailenter_lik = lik;
  const double tailenter_homr = homr;
  const double tailenter_homc = homc;
  const double tailenter_hets = hets;
  while (lik <= one_minus_scaled_eps) {
    tail_sum += lik;
    homr += 1;
    homc += 1;
    lik *= hets * (hets - 1) / (4 * homr * homc);
    hets -= 2;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lik < 2 - one_minus_scaled_eps) {
    const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
    const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &lik);
    if (cmp_result <= 0) {
      tail_sum += lik;
      if (midp && (cmp_result == 0)) {
        tail_sum -= 0.5;
      }
    }
  }
  // Sum away from center, until sums stop changing.
  lik = tailenter_lik;
  homr = tailenter_homr;
  homc = tailenter_homc;
  hets = tailenter_hets;
  while (1) {
    hets += 2;
    lik *= 4 * homr * homc / (hets * (hets - 1));
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
    homr -= 1;
    homc -= 1;
  }
  return starting_lnprob + log(tail_sum);
}

// 2^{-83} bias to give plink 1.9-style exact tests maximum ability to
// determine tiny p-values.  (~2^{-53} is necessary to take advantage of
// denormalized small numbers, then allow tail sum to be up to 2^30.  ...okay,
// HweThresh[Midp]() is not responsible for denormal values of thresh, and
// plink2 now just flushes denormals to zero.  But configuring this
// constant to be compatible with them doesn't cost us anything.)
static const double kExactTestBias = k2m50 / (1LL << 33);

uint32_t HweThresh(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh) {
  // Threshold-test-only version of HweLnP() which is usually able to exit
  // from the calculation earlier.  Assumes DBL_MIN <= pval_thresh <= 1 (note
  // that some older versions of this function didn't handle pval_thresh=1
  // correctly, and plink2 still avoids calling this with pval_thresh=1).
  // Returns 0 if these counts are close enough to Hardy-Weinberg equilibrium,
  // 1 otherwise.
  //
  // Suppose, for definiteness, that the number of observed hets is no less
  // than expectation.  (Same ideas apply for the other case.)  We proceed as
  // follows:
  // - Sum the *relative* likelihoods of more-likely smaller het counts.
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
    return 0;
  }
  const int32_t rare_copies = 2 * obs_homr + obs_hets;
  double hets_t2 = obs_hets;  // tail 2
  double homr_t2 = obs_homr;
  double homc_t2 = obs_homc;

  double tail_sum1 = kExactTestBias;
  double center_sum = 0;
  double lik2 = tail_sum1;
  double tail_sum2 = 0;

  // const double center_div_tail_thresh = (1 - pval_thresh) / pval_thresh;
  // An initial upper bound on the tail sum is useful, since it lets us
  // report test failure before summing the entire center.  We use the
  // trivial bound of 1 + floor(rare_copies / 2): that's the total number
  // of possible het counts, and the relative likelihood for each count must be
  // <= 1 if it's in the tail.
  const double center_div_tail_thresh = (1 - pval_thresh) / pval_thresh;
  const double center_sum_exit_thresh = u31tod(1 + (rare_copies >> 1)) * (center_div_tail_thresh * kExactTestBias);
  double scaled_one_plus_eps = kExactTestBias * (1 + k2m52);

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
      return 0;
    }

    // het_probs[hets] = 1
    // het_probs[hets - 2] = het_probs[hets] * hets * (hets - 1) / (4 * (homr + 1) * (homc + 1))
    do {
      homr_t2 += 1;
      homc_t2 += 1;
      lik2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      hets_t2 -= 2;
      scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
      if (lik2 < scaled_one_plus_eps) {
        if (lik2 <= 2 * kExactTestBias - scaled_one_plus_eps) {
          tail_sum2 = lik2;
          break;
        }
        double unshifted_lik = lik2 * (1.0 / kExactTestBias);
        const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
        const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, nullptr, &unshifted_lik);
        lik2 = unshifted_lik * kExactTestBias;
        if (cmp_result <= 0) {
          tail_sum2 = lik2;
          break;
        }
        // HweCompare() could have computed unshifted_like using
        // exp(lnprobv_diff).  We are conservatively assuming that exp()
        // introduces up to 1.5 ULP error, and lnprobv_diff may be off by
        // slightly over 0.5 ULP.
        scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
      }
      center_sum += lik2;
      if (center_sum >= center_sum_exit_thresh) {
        return 1;
      }
    } while (hets_t2 > 1);
    // hets_t2 guaranteed to be nonnegative on loop exit, so ratio
    // calculation works

    // This is NaN when pval_thresh=1, so we write the next if-condition to be
    // true on NaN.
    const double tail_sum_exit_thresh = center_sum / center_div_tail_thresh;
    if (!(tail_sum1 + tail_sum2 < tail_sum_exit_thresh)) {
      return 0;
    }
    // c + cr + cr^2 + ... = c/(1-r), which is an upper bound for the tail sum
    const double ratio = (hets_t2 * (hets_t2 - 1)) / (4 * (homr_t2 + 1) * (homc_t2 + 1));
    const double tail2_ceil = tail_sum2 / (1 - ratio);
    double hets_t1 = obs_hets + 2;
    double homr_t1 = obs_homr;
    double homc_t1 = obs_homc;
    // ratio for the other tail
    double lik1 = (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
    const double tail1_ceil = tail_sum1 / (1 - lik1);
    if (tail1_ceil + tail2_ceil < tail_sum_exit_thresh) {
      return 1;
    }
    lik1 *= tail_sum1;
    tail_sum1 += lik1;

    if (obs_homr > 1) {
      // het_probs[hets + 2] = het_probs[hets] * 4 * homr * homc / ((hets + 2) * (hets + 1))
      const double tail_sum1_exit_thresh = tail_sum_exit_thresh - tail_sum2;
      while (1) {
        hets_t1 += 2;
        homr_t1 -= 1;
        homc_t1 -= 1;
        lik1 *= (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
        const double preadd = tail_sum1;
        tail_sum1 += lik1;
        if (tail_sum1 >= tail_sum1_exit_thresh) {
          return 0;
        }
        // homr_t1 == 1 check isn't necessary for correctness, but it provides
        // a noticeable speedup in my testing on real data.
        if ((tail_sum1 == preadd) || (homr_t1 == 1)) {
          break;
        }
      }
    }
    if (tail_sum1 + tail2_ceil < tail_sum_exit_thresh) {
      return 1;
    }
    const double tail_sum2_exit_thresh = tail_sum_exit_thresh - tail_sum1;
    while (1) {
      homr_t2 += 1;
      homc_t2 += 1;
      lik2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      const double preadd = tail_sum2;
      tail_sum2 += lik2;
      if (tail_sum2 >= tail_sum2_exit_thresh) {
        return 0;
      }
      if (tail_sum2 == preadd) {
        return 1;
      }
      hets_t2 -= 2;
    }
  }
  // tail 1 = lower hets
  if (!obs_homr) {
    return 0;
  }
  do {
    hets_t2 += 2;
    lik2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
    if (lik2 < scaled_one_plus_eps) {
      if (lik2 <= 2 * kExactTestBias - scaled_one_plus_eps) {
        tail_sum2 = lik2;
        break;
      }
      double unshifted_lik = lik2 * (1.0 / kExactTestBias);
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
      const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, nullptr, &unshifted_lik);
      lik2 = unshifted_lik * kExactTestBias;
      if (cmp_result <= 0) {
        tail_sum2 = lik2;
        break;
      }
      scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
    }
    center_sum += lik2;
    if (center_sum >= center_sum_exit_thresh) {
      return 1;
    }
  } while (homr_t2 > 0);
  // homr_t2 guaranteed to be nonnegative on loop exit, so ratio
  // calculation works

  const double tail_sum_exit_thresh = center_sum / center_div_tail_thresh;
  if (!(tail_sum1 + tail_sum2 < tail_sum_exit_thresh)) {
    return 0;
  }
  const double ratio = (4 * homr_t2 * homc_t2) / ((hets_t2 + 2) * (hets_t2 + 1));
  const double tail2_ceil = tail_sum2 / (1 - ratio);
  double hets_t1 = obs_hets;
  double homr_t1 = obs_homr + 1;
  double homc_t1 = obs_homc + 1;
  double lik1 = (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
  const double tail1_ceil = tail_sum1 / (1 - lik1);
  if (tail1_ceil + tail2_ceil < tail_sum_exit_thresh) {
    return 1;
  }
  lik1 *= tail_sum1;
  tail_sum1 += lik1;

  if (obs_hets >= 4) {
    const double tail_sum1_exit_thresh = tail_sum_exit_thresh - tail_sum2;
    while (1) {
      hets_t1 -= 2;
      homr_t1 += 1;
      homc_t1 += 1;
      lik1 *= (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
      const double preadd = tail_sum1;
      tail_sum1 += lik1;
      if (tail_sum1 >= tail_sum1_exit_thresh) {
        return 0;
      }
      if ((tail_sum1 == preadd) || (hets_t1 < 4)) {
        break;
      }
    }
  }
  if (tail_sum1 + tail2_ceil < tail_sum_exit_thresh) {
    return 1;
  }
  const double tail_sum2_exit_thresh = tail_sum_exit_thresh - tail_sum1;
  while (1) {
    hets_t2 += 2;
    lik2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    const double preadd = tail_sum2;
    tail_sum2 += lik2;
    if (tail_sum2 >= tail_sum2_exit_thresh) {
      return 0;
    }
    if (tail_sum2 == preadd) {
      return 1;
    }
  }
}

uint32_t HweThreshMidp(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, double pval_thresh) {
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
    return 0;
  }
  int32_t rare_copies = 2 * obs_homr + obs_hets;
  double hets_t2 = obs_hets;  // tail 2
  double homr_t2 = obs_homr;
  double homc_t2 = obs_homc;
  double tail_sum1 = kExactTestBias * 0.5;
  double center_sum = tail_sum1;
  double lik2 = kExactTestBias;
  double tail_sum2 = 0;
  const double center_div_tail_thresh = (1 - pval_thresh) / pval_thresh;
  const double center_sum_exit_thresh = u31tod(1 + (rare_copies >> 1)) * (center_div_tail_thresh * kExactTestBias);
  double scaled_one_plus_eps = kExactTestBias * (1 + k2m52);
  if (obs_hets * genotypes2 > rare_copies * (genotypes2 - rare_copies)) {
    if (obs_hets < 2) {
      return 0;
    }
    do {
      homr_t2 += 1;
      homc_t2 += 1;
      lik2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      hets_t2 -= 2;
      scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
      if (lik2 < scaled_one_plus_eps) {
        if (lik2 <= 2 * kExactTestBias - scaled_one_plus_eps) {
          tail_sum2 = lik2;
          break;
        }
        double lik = lik2 * (1.0 / kExactTestBias);
        const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
        const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, nullptr, &lik);
        lik2 = lik * kExactTestBias;
        if (cmp_result <= 0) {
          if (cmp_result == 0) {
            tail_sum2 = tail_sum1;
            center_sum += tail_sum1;
          } else {
            tail_sum2 = lik2;
          }
          break;
        }
        scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
      }
      center_sum += lik2;
      if (center_sum >= center_sum_exit_thresh) {
        return 1;
      }
    } while (hets_t2 > 1);

    const double tail_sum_exit_thresh = center_sum / center_div_tail_thresh;
    if (!(tail_sum1 + tail_sum2 < tail_sum_exit_thresh)) {
      return 0;
    }
    const double ratio = (hets_t2 * (hets_t2 - 1)) / (4 * (homr_t2 + 1) * (homc_t2 + 1));
    // this needs to work in both the tie and no-tie cases
    const double tail2_ceil = prefer_fma(lik2, ratio / (1 - ratio), tail_sum2);
    double hets_t1 = obs_hets + 2;
    double homr_t1 = obs_homr;
    double homc_t1 = obs_homc;
    double lik1 = (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
    // always a tie here
    const double tail1_ceil = (tail_sum1 * 2) / (1 - lik1) - tail_sum1;
    if (tail1_ceil + tail2_ceil < tail_sum_exit_thresh) {
      return 1;
    }
    lik1 *= tail_sum1 * 2;
    tail_sum1 += lik1;

    if (obs_homr > 1) {
      const double tail_sum1_exit_thresh = tail_sum_exit_thresh - tail_sum2;
      while (1) {
        hets_t1 += 2;
        homr_t1 -= 1;
        homc_t1 -= 1;
        lik1 *= (4 * homr_t1 * homc_t1) / (hets_t1 * (hets_t1 - 1));
        const double preadd = tail_sum1;
        tail_sum1 += lik1;
        if (tail_sum1 >= tail_sum1_exit_thresh) {
          return 0;
        }
        if ((tail_sum1 == preadd) || (homr_t1 == 1)) {
          break;
        }
      }
    }
    if (tail_sum1 + tail2_ceil < tail_sum_exit_thresh) {
      return 1;
    }
    const double tail_sum2_exit_thresh = tail_sum_exit_thresh - tail_sum1;
    while (1) {
      homr_t2 += 1;
      homc_t2 += 1;
      lik2 *= (hets_t2 * (hets_t2 - 1)) / (4 * homr_t2 * homc_t2);
      const double preadd = tail_sum2;
      tail_sum2 += lik2;
      if (tail_sum2 >= tail_sum2_exit_thresh) {
        return 0;
      }
      if (tail_sum2 == preadd) {
        return 1;
      }
      hets_t2 -= 2;
    }
  }
  if (!obs_homr) {
    return 0;
  }
  do {
    hets_t2 += 2;
    lik2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    scaled_one_plus_eps += kExactTestBias * 2 * k2m52;
    if (lik2 < scaled_one_plus_eps) {
      if (lik2 <= 2 - scaled_one_plus_eps) {
        tail_sum2 = lik2;
        break;
      }
      double lik = lik2 * (1.0 / kExactTestBias);
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr_t2);
      const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, nullptr, &lik);
      lik2 = lik * kExactTestBias;
      if (cmp_result <= 0) {
        if (cmp_result == 0) {
          tail_sum2 = tail_sum1;
          center_sum += tail_sum1;
        } else {
          tail_sum2 = lik2;
        }
        break;
      }
      scaled_one_plus_eps = kExactTestBias * (1 + 3 * k2m52);
    }
    center_sum += lik2;
    if (center_sum >= center_sum_exit_thresh) {
      return 1;
    }
  } while (homr_t2 > 0);

  const double tail_sum_exit_thresh = center_sum / center_div_tail_thresh;
  if (!(tail_sum1 + tail_sum2 < tail_sum_exit_thresh)) {
    return 0;
  }
  const double ratio = (4 * homr_t2 * homc_t2) / ((hets_t2 + 2) * (hets_t2 + 1));
  const double tail2_ceil = prefer_fma(lik2, ratio / (1 - ratio), tail_sum2);
  double hets_t1 = obs_hets;
  double homr_t1 = obs_homr + 1;
  double homc_t1 = obs_homc + 1;
  double lik1 = (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
  const double tail1_ceil = (2 * tail_sum1) / (1 - lik1) - tail_sum1;
  lik1 *= 2 * tail_sum1;
  tail_sum1 += lik1;

  if (tail1_ceil + tail2_ceil < tail_sum_exit_thresh) {
    return 1;
  }
  if (obs_hets >= 4) {
    const double tail_sum1_exit_thresh = tail_sum_exit_thresh - tail_sum2;
    while (1) {
      hets_t1 -= 2;
      homr_t1 += 1;
      homc_t1 += 1;
      lik1 *= (hets_t1 * (hets_t1 - 1)) / (4 * homr_t1 * homc_t1);
      const double preadd = tail_sum1;
      tail_sum1 += lik1;
      if (tail_sum1 >= tail_sum1_exit_thresh) {
        return 0;
      }
      if ((tail_sum1 == preadd) || (hets_t1 < 4)) {
        break;
      }
    }
  }
  if (tail_sum1 + tail2_ceil < tail_sum_exit_thresh) {
    return 1;
  }
  const double tail_sum2_exit_thresh = tail_sum_exit_thresh - tail_sum1;
  while (1) {
    hets_t2 += 2;
    lik2 *= (4 * homr_t2 * homc_t2) / (hets_t2 * (hets_t2 - 1));
    homr_t2 -= 1;
    homc_t2 -= 1;
    const double preadd = tail_sum2;
    tail_sum2 += lik2;
    if (tail_sum2 >= tail_sum2_exit_thresh) {
      return 0;
    }
    if (tail_sum2 == preadd) {
      return 1;
    }
  }
}

uint32_t HweThreshLnMain(int32_t obs_hets, int32_t obs_hom1, int32_t obs_hom2, int32_t midp, double ln_thresh) {
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
    return 0;
  }

  // 1. Compute log-probability of starting table.  This may be high enough on
  //    its own to immediately return 0.
  //    If probability is lower than threshold / <total # of tables>, we can
  //    immediately return 1.
  // 2. Determine tailsum we must hit to return 0.
  // 3. The rest follows HweLnP(), except with an extra geometric-series-based
  //    early-exit attempt near the end.
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_add3_lfacts(rare_ct, sample_ct, allele_ctd - rare_ct),
            ddr_lfact(allele_ctd));
  double homr = obs_homr;
  double homc = obs_homc;
  dd_real starting_lnprobv_ddr =
    ddr_sub(ddr_muld(_ddr_log2, hets),
            ddr_add3_lfacts(homr, hets, homc));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  if (ln_thresh <= starting_lnprob - midp * kLn2) {
    return 0;
  }
  const double max_homr = S_CAST(double, rare_ct >> 1);
  if (ln_thresh > starting_lnprob + log(max_homr + 1 - midp * 0.5)) {
    return 1;
  }

  const double c_minus_r = homc - homr;
  // This should be in (0.5, 2^30].
  const double tail_thresh = exp(ln_thresh - starting_lnprob);
  double tail_sum = 1 - midp * 0.5;
  double lik = 1;
  if (hets > cmodal_nhet) {
    // No center-sum (or p=1) code path, since it doesn't make sense to choose
    // a HWE threshold that makes these relevant; we should have already
    // returned 0.
    // (So we can't assume obs_hets >= 2.)
    while (1) {
      hets += 2;
      lik *= (4 * homr * homc) / (hets * (hets - 1));
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        break;
      }
      if (tail_sum >= tail_thresh) {
        return 0;
      }
      homr -= 1;
      homc -= 1;
    }
    {
      const double delta = 0.5 * (hets + obs_hets) - cmodal_nhet;
      homr = 0.5 * (homr + obs_homr) + delta;
      homr = ceil_limit(homr, max_homr);
    }
    while (1) {
      hets = rare_ct - homr * 2;
      homc = homr + c_minus_r;
      const dd_real lnprobv_ddr =
        ddr_sub(ddr_muld(_ddr_log2, hets),
                ddr_add3_lfacts(homr, hets, homc));
      const double lnprobv_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
      if (lnprobv_diff >= k2m60) {
        if (homr == max_homr) {
          return 1;
        }
        const double lnprobv_deriv = log(hets * (hets - 1) / (4 * (homr + 1) * (homc + 1)));
        homr += ceil(-lnprobv_diff / lnprobv_deriv);
        if (homr > max_homr) {
          homr = max_homr;
        }
      } else if (lnprobv_diff > -62 * kLn2) {
        lik = exp(lnprobv_diff);
        break;
      } else {
        const double lnprobv_deriv = log((hets + 2) * (hets + 1) / (4 * homr * homc));
        homr -= S_CAST(int64_t, lnprobv_diff / lnprobv_deriv);
        assert(homr >= 0);
      }
    }
    double one_minus_scaled_eps = 1 - 3 * k2m52;
    const double tailenter_lik = lik;
    const double tailenter_homr = homr;
    const double tailenter_homc = homc;
    const double tailenter_hets = hets;
    while (lik <= one_minus_scaled_eps) {
      tail_sum += lik;
      if (tail_sum >= tail_thresh) {
        return 0;
      }
      hets += 2;
      lik *= (4 * homr * homc) / (hets * (hets - 1));
      homr -= 1;
      homc -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lik < 2 - one_minus_scaled_eps) {
      const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
      const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &lik);
      if (cmp_result <= 0) {
        if (cmp_result == 0) {
          tail_sum += 1 - midp * 0.5;
        } else {
          tail_sum += lik;
        }
        if (tail_sum >= tail_thresh) {
          return 0;
        }
      }
    }
    // We're down to one tail that can be tightly bounded by a geometric
    // series.  (ratio is always decreasing)
    // c + cr + cr^2 + ... = c/(1-r)
    lik = tailenter_lik;
    homr = tailenter_homr;
    homc = tailenter_homc;
    hets = tailenter_hets;

    homr += 1;
    homc += 1;
    const double cur_ratio = (hets * (hets - 1)) / (4 * homr * homc);
    lik *= cur_ratio;
    const double remaining_ceil = tailenter_lik / (1 - cur_ratio);
    if (tail_sum + remaining_ceil < tail_thresh) {
      return 1;
    }
    while (1) {
      const double preadd = tail_sum;
      tail_sum += lik;
      if (tail_sum == preadd) {
        return 1;
      }
      if (tail_sum >= tail_thresh) {
        return 0;
      }
      hets -= 2;
      homr += 1;
      homc += 1;
      lik *= (hets * (hets - 1)) / (4 * homr * homc);
    }
    // unreachable
  }
  while (1) {
    homr += 1;
    homc += 1;
    lik *= hets * (hets - 1) / (4 * homr * homc);
    hets -= 2;
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      break;
    }
    if (tail_sum >= tail_thresh) {
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
              ddr_add3_lfacts(homr, hets, homc));
    const double lnprobv_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    if (lnprobv_diff >= k2m60) {
      if (homr <= 0) {
        return 1;
      }
      const double lnprobv_deriv = log(4 * homr * homc / ((hets + 2) * (hets + 1)));
      homr -= ceil(-lnprobv_diff / lnprobv_deriv);
      if (homr < 0) {
        homr = 0;
      }
    } else if (lnprobv_diff > -62 * kLn2) {
      lik = exp(lnprobv_diff);
      break;
    } else {
      const double lnprobv_deriv = log(4 * (homr + 1) * (homc + 1) / (hets * (hets - 1)));
      homr += S_CAST(int64_t, lnprobv_diff / lnprobv_deriv);
      assert(homr <= max_homr);
    }
  }
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double tailenter_lik = lik;
  const double tailenter_homr = homr;
  const double tailenter_homc = homc;
  const double tailenter_hets = hets;
  while (lik <= one_minus_scaled_eps) {
    tail_sum += lik;
    if (tail_sum >= tail_thresh) {
      return 0;
    }
    homr += 1;
    homc += 1;
    lik *= hets * (hets - 1) / (4 * homr * homc);
    hets -= 2;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lik < 2 - one_minus_scaled_eps) {
    const int32_t hom_decr = obs_homr - S_CAST(int32_t, homr);
    const intptr_t cmp_result = HweCompare(obs_hets, obs_homr, obs_homc, hom_decr, &starting_lnprobv_ddr, &lik);
    if (cmp_result <= 0) {
      if (cmp_result == 0) {
        tail_sum += 1 - midp * 0.5;
      } else {
        tail_sum += lik;
      }
      if (tail_sum >= tail_thresh) {
        return 0;
      }
    }
  }
  lik = tailenter_lik;
  homr = tailenter_homr;
  homc = tailenter_homc;
  hets = tailenter_hets;

  hets += 2;
  const double cur_ratio = 4 * homr * homc / (hets * (hets - 1));
  lik *= cur_ratio;
  const double remaining_ceil = lik / (1 - cur_ratio);
  if (tail_sum + remaining_ceil < tail_thresh) {
    return 1;
  }
  while (1) {
    const double preadd = tail_sum;
    tail_sum += lik;
    if (tail_sum == preadd) {
      return 1;
    }
    if (tail_sum >= tail_thresh) {
      return 0;
    }
    homr -= 1;
    homc -= 1;
    hets += 2;
    lik *= 4 * homr * homc / (hets * (hets - 1));
  }
}

#ifdef __cplusplus
}
#endif
