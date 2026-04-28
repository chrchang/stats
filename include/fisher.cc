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

#include "fisher.h"

#include <assert.h>
#include <math.h>

#include "plink2_float.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// *cmp_resultp is set to positive value if m22 = obs_m22 + m22_incr has higher
// likelihood than m22 = obs_m22, 0 if identical likelihood, and negative value
// if lower likelihood.
// Error is returned iff malloc fails.
// If neg_numer_ddr has not been computed yet, set its x[0] to DBL_MAX; it will
// be filled in if necessary.
BoolErr Fisher22Compare(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m21, uint32_t obs_m22, int32_t m22_incr, dd_real* neg_numer_ddr_ptr, intptr_t* cmp_resultp, double* dbl_ptr) {
  // Fisher 2x2 likelihood is
  //
  //   (m11+m12)! (m21+m22)! (m11+m21)! (m12+m22)!
  //   -------------------------------------------
  //     m11! m12! m21! m22! (m11+m12+m21+m22)!
  //
  // so likelihood ratio of interest is
  //
  //           obs_m11! obs_m12! obs_m21! obs_m22!
  //   ---------------------------------------------------
  //   (obs_m11+j)! (obs_m12-j)! (obs_m21-j)! (obs_m22+j)!
  //
  // where j=m22_incr.
  //
  // Note that HWE kind of maps to this, via
  //   m11 := obs_hets*0.5
  //   m12 := obs_hom1
  //   m21 := obs_hom2
  //   m22 := (obs_hets-1)*0.5
  uint32_t numer_factorial_args[4];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m21;
  numer_factorial_args[3] = obs_m22;
  uint32_t denom_factorial_args[4];
  denom_factorial_args[0] = obs_m11 + m22_incr;
  denom_factorial_args[1] = obs_m12 - m22_incr;
  denom_factorial_args[2] = obs_m21 - m22_incr;
  denom_factorial_args[3] = obs_m22 + m22_incr;

  mp_limb_t* gmp_wkspace = nullptr;
  uintptr_t gmp_wkspace_limb_ct = 0;
  BoolErr reterr = CompareFactorialProducts(4, 0, 0, numer_factorial_args, denom_factorial_args, neg_numer_ddr_ptr, &gmp_wkspace, &gmp_wkspace_limb_ct, cmp_resultp, dbl_ptr);
  free_cond(gmp_wkspace);
  return reterr;
}

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^31.
// They're defined as int32_ts instead of uint32_ts since signed int <->
// floating-point conversions are sometimes faster than the same-width unsigned
// int <-> floating-point conversions.  (Yes, uint32 -> float64 should be fine
// on 64-bit systems, but we may as well write this to also be efficient on
// 32-bit systems when there's no real drawback to doing so.)
BoolErr Fisher22LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, int32_t midp, double* resultp) {
  // Normalize: m11 >= m22, m12 >= m21, m11*m22 < m12*m21.
  // Note that the first two are reversed from PLINK 1.9, to get rid of
  // spurious index differences between Fisher22 and Fisher23.
  if (obs_m11 < obs_m22) {
    swap_i32(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i32(&obs_m12, &obs_m21);
  }
  if (S_CAST(int64_t, obs_m11) * obs_m22 > S_CAST(int64_t, obs_m12) * obs_m21) {
    swap_i32(&obs_m11, &obs_m12);
    swap_i32(&obs_m21, &obs_m22);
  }
  if (!midp) {
    // Fast path for p=1.
    if ((obs_m11 + 1LL) * (obs_m22 + 1) == S_CAST(int64_t, obs_m12) * obs_m21) {
      *resultp = 0;
      return 0;
    }
  }
  // Iterate outward to floating-point precision limit.

  // I experimented with making more int32 -> double casts explicit, but
  // concluded that it was hurting readability more than it was helping.
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  double lastp = 1;
  double tailp = 1 - 0.5 * midp;
  while (1) {
    m12 += 1;
    m21 += 1;
    lastp *= (m11 * m22) / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
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
  // As with HweLnP(), we note that if we're within 172 steps of the mode and
  // the starting relative-likelihood is normalized to 1, the modal
  // relative-likelihood can be loosely bounded above by
  //   ((172^172) / 172!)^4 ~= 5.3e+292
  // which leaves enough headroom to accumulate the rest of the center-sum and
  // represent intermediate values without overflowing.
  if ((obs_m11 + 172LL) * (obs_m22 + 172LL) >= (obs_m12 - 172LL) * (obs_m21 - 172LL)) {
    lastp = 1;
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    dd_real starting_lnprobv_ddr = {{DBL_MAX, 0.0}};
    double centerp = 0.5 * midp;
    while (1) {
      m11 += 1;
      m22 += 1;
      lastp *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      // Number of center contingency tables is maximized with obs_m22 = 0,
      // modal_m22 = 172, other values large.
      // Since 1 + 1/2 + ... + 1/172 < 1/173 + ... + 1/53000, we're limited to
      // ~53000 tables.  Each lastp update involves 4 operations which can each
      // introduce up to 0.5 ULP relative error under the default rounding
      // mode.
      if (lastp < 1 + 53000 * 2 * k2m52) {
        if (lastp <= 1 - 53000 * 2 * k2m52) {
          tailp += lastp;
          break;
        }
        // Near-tie.  True value of lastp can be greater than, equal to, or
        // less than 1.
        const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
        intptr_t cmp_result;
        if (unlikely(Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
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
      m11 += 1;
      m22 += 1;
      lastp *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
    }
    *resultp = log(tailp / (tailp + centerp));
    return 0;
  }
  dd_real starting_lnprobv_ddr =
    ddr_negate(ddr_add4_lfacts(obs_m11, obs_m12, obs_m21, obs_m22));

  // Now we want to jump near the other tail, without evaluating that many
  // contingency table log-likelihoods along the way.
  //
  // Each full log-likelihood evaluation requires 4 ddr_lfact() calls.  Since
  // they are now performed with extra precision, they require hundreds of
  // floating-point operations, so we want to limit ourselves to 1-2 full
  // evaluations most of the time.  (Possible todo: use lower-accuracy Lfact()
  // to jump around, followed by ddr_lfact() when exiting the loop.  Should be
  // an easy performance win, but there's a complexity cost so I'll wait until
  // I see a scenario where this branch executes frequently...)
  //
  // The current heuristic starts by reflecting (obs_m21 + m21) * 0.5 across
  // the mode, performing a full log-likelihood check at an adjacent valid
  // point.  (It is convenient to focus on m21 here, since m21=0 corresponds to
  // the outermost table on this tail.)  Hopefully we find that we're in
  // (starting_lnprob - 62 * kLn2, starting_lnprob], so we're at or near a
  // table that actually contributes to the tail-sum.  (This window is chosen
  // to be wide enough to guarantee that at least one point falls inside when
  // obs_m11 + obs_m12 + obs_m21 + obs_m22 < 2^31.)
  //
  // If not, we jump again, using Newton's method.
  // If m21 is too high (i.e. current log-likelihood is too high), decreasing
  // m12 by 1 would multiply the likelihood by
  //   (m11 + 1) * (m22 + 1) / (m12 * m21)
  // If m12 is too low, increasing m12 by 1 would multiply the likelihood by
  //   (m12 + 1) * (m21 + 1) / (m11 * m22)
  // We use the negative-log of the first expression as the Newton's method
  // f'(x) when we're jumping to lower m12, and the log of the second
  // expression when we're jumping to higher m12.
  // f''(x) is always negative, so we can aim for starting_lnprob instead of
  // the middle of the interval.

  const double m12_minus_m21 = obs_m12 - obs_m21;
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx1 = obs_m11 + obs_m21;
  const double mxx = m1x + m2x;
  {
    // x=modal_m21 satisfies
    //    (m2x - x) * (mx1 - x) = x * (x + m12_minus_m21)
    // -> (x - m2x) * (x - mx1) = x * (x + m12_minus_m21)
    // -> x^2 + x*(-m2x - mx1) + m2x*mx1 = x^2 + x*m12_minus_m21
    // -> m2x*mx1 = x*(m12_minus_m21 + m2x + mx1)
    // -> x = m2x*mx1 / mxx
    const double modal_m21 = m2x * mx1 / mxx;
    m21 = 2 * modal_m21 - (m21 + obs_m21) * 0.5;
    // Round down (to guarantee we've actually moved to the other side of the
    // mode) and clamp.
    m21 = S_CAST(int32_t, m21);
    if (m21 < 0) {
      m21 = 0;
    }
  }
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_add4_lfacts(m1x, m2x, mx1, obs_m12 + obs_m22),
            ddr_lfact(mxx));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  while (1) {
    m11 = mx1 - m21;
    m12 = m12_minus_m21 + m21;
    m22 = m2x - m21;
    const dd_real lnprobv_ddr =
      ddr_negate(ddr_add4_lfacts(m11, m12, m21, m22));
    const double lnprob_diff = ddr_sub(lnprobv_ddr, starting_lnprobv_ddr).x[0];
    // Could tighten this threshold further; I haven't performed a careful
    // error analysis yet but CompareFactorialProducts() includes a plausible
    // assumption that 2^{-60} is safe.  But code is correct as long as we're
    // guaranteed to enter the "lastp < 2 - one_minus_scaled_eps" branch for
    // positive lnprob_diff.
    if (lnprob_diff >= k2m53) {
      if (m21 == 0) {
        // All tables on this tail have higher likelihood than the starting
        // table.  Exit.
        if (midp) {
          tailp -= 0.5;
        }
        *resultp = starting_lnprob + log(tailp);
        return 0;
      }
      const double ll_deriv = log(m12 * m21 / ((m11 + 1) * (m22 + 1)));
      // Round up, to guarantee that we make progress.
      // (lnprob_diff is positive and ll_deriv is negative.)
      // This may overshoot.  But the function is guaranteed to terminate
      // because we never overshoot (and we do always make progress on each
      // step) once we're on the other side.
      m21 -= ceil_smalleps(-lnprob_diff / ll_deriv);
      if (m21 < 0) {
        m21 = 0;
      }
    } else if (lnprob_diff > -62 * kLn2) {
      lastp = exp(lnprob_diff);
      break;
    } else {
      const double ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
      // Round down, to guarantee we don't overshoot.
      // We're guaranteed to make progress, since lnprob_diff <= -62 * log(2),
      // m11 * m22 < 2^62, and (m12 + 1) * (m21 + 1) >= 1.
      m21 += S_CAST(int64_t, lnprob_diff / ll_deriv);
    }
  }
  // Sum toward center, until lastp >= 1.
  // lastp should be accurate to 3 ULP as we enter this loop (max 1.5 ULP
  // observed error from exp, tiny bit over 0.5 from lnprob_diff), so near-tie
  // detection can use a tight epsilon here.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  const double lastp_tail = lastp;
  const double m11_tail = m11;
  const double m12_tail = m12;
  const double m21_tail = m21;
  const double m22_tail = m22;
  while (lastp <= one_minus_scaled_eps) {
    tailp += lastp;
    m12 += 1;
    m21 += 1;
    lastp *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lastp < 2 - one_minus_scaled_eps) {
    const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
    intptr_t cmp_result;
    if (unlikely(Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, &starting_lnprobv_ddr, &cmp_result, &lastp))) {
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
  m11 = m11_tail;
  m12 = m12_tail;
  m21 = m21_tail;
  m22 = m22_tail;
  while (1) {
    m11 += 1;
    m22 += 1;
    lastp *= m12 * m21 / (m11 * m22);
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
    m12 -= 1;
    m21 -= 1;
  }
  *resultp = starting_lnprob + log(tailp);
  return 0;
}

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^31.
// Just returns 0 if the log-p-value > 1 - 2^{-54} since additional precision
// in that direction is expected to be irrelevant.  (Reverse the test direction
// when you do actually want that precision.)
double Fisher22OneSidedLnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, uint32_t m11_is_greater_alt, int32_t midp) {
  // Normalize.
  if (obs_m11 < obs_m22) {
    swap_i32(&obs_m11, &obs_m22);
  }
  if (obs_m12 < obs_m21) {
    swap_i32(&obs_m12, &obs_m21);
  }
  // Flipping m11<->m12 and m21<->m22 also flips the direction of the
  // alternative hypothesis.  So we flip on m11-is-greater alternative
  // hypothesis here to allow the rest of the code to assume m11-is-less.
  if (m11_is_greater_alt) {
    swap_i32(&obs_m11, &obs_m12);
    swap_i32(&obs_m21, &obs_m22);
  }
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  if (S_CAST(int64_t, obs_m11) * obs_m22 >= S_CAST(int64_t, obs_m12) * obs_m21) {
    // We're at or to the right of the mode.
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know the p-value > 1 - 2^{-54} (at which point
    // we return 0), or remaining left likelihoods are smaller than the
    // precision limit.
    double right_sum = m12 * m21 / ((m11 + 1) * (m22 + 1));
    // r + r^2 + ... = r / (1-r)
    // const double right_upper_bound = right_sum / (1 - right_sum);

    // Rescale our starting lastp so that we overflow to INFINITY when we'd
    // want to early-exit and return 0; this saves us a comparison in the loop.
    const double left_rescale = (DBL_MAX / (1LL << 54)) * (1 - right_sum) / right_sum;
    double lastp = left_rescale;
    double left_sum = left_rescale;
    while (1) {
      m12 += 1;
      m21 += 1;
      lastp *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      const double preaddp = left_sum;
      left_sum += lastp;
      if (left_sum == preaddp) {
        break;
      }
    }
    if (left_sum == INFINITY) {
      return 0;
    }
    left_sum /= left_rescale;

    // Now compute the right-sum to the precision limit.
    m11 = obs_m11 + 1.0;
    m12 = obs_m12 - 1;
    m21 = obs_m21 - 1;
    m22 = obs_m22 + 1;
    lastp = right_sum;
    while (1) {
      m11 += 1;
      m22 += 1;
      lastp *= m12 * m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preaddp = right_sum;
      right_sum += lastp;
      if (right_sum == preaddp) {
        break;
      }
    }
    return log1p((-0.5 * midp - right_sum) / (right_sum + left_sum));
  }
  // We're to the left of the mode, and are responsible for tiny p-values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise... if the problem instance isn't *that* large, we could use
  // Lfact() to compute the starting log-likelihood and eat a big catastrophic
  // cancellation error, but I'll start with just the slow-and-accurate
  // calculation.
  const double m11_minus_m22 = obs_m11 - obs_m22;
  const double m1x = obs_m11 + obs_m12;
  const double m2x = obs_m21 + obs_m22;
  const double mx2 = obs_m12 + obs_m22;
  const double mxx = m1x + m2x;
  const double modal_m22 = m2x * mx2 / mxx;
  if (modal_m22 <= S_CAST(double, 172LL + obs_m22)) {
    double lastp = 1;
    double right_sum = 0;
    while (1) {
      m11 += 1;
      m22 += 1;
      lastp *= m12 * m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preaddp = right_sum;
      right_sum += lastp;
      if (right_sum == preaddp) {
        break;
      }
    }
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    lastp = 1;
    double left_sum = 1;
    while (1) {
      m12 += 1;
      m21 += 1;
      lastp *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      const double preaddp = left_sum;
      left_sum += lastp;
      if (left_sum == preaddp) {
        break;
      }
    }
    return log((left_sum - 0.5 * midp) / (left_sum + right_sum));
  }
  const dd_real starting_lnprobv_ddr =
    ddr_negate(ddr_add4_lfacts(obs_m11, obs_m12, obs_m21, obs_m22));
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_add4_lfacts(m1x, m2x, obs_m11 + obs_m21, mx2),
            ddr_lfact(mxx));
  const double starting_lnprob = ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  double lastp = 1;
  double left_sum = 1;
  while (1) {
    m12 += 1;
    m21 += 1;
    lastp *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preaddp = left_sum;
    left_sum += lastp;
    if (left_sum == preaddp) {
      break;
    }
  }
  return starting_lnprob + log(left_sum - 0.5 * midp);
}

// Switch between log- and regular representations at kSwitchThresh.
// 2^890 leaves enough headroom for at least 4 more multiplies by obs_total.
// (Useful to reduce this to e.g. 2^150 when testing correctness, though it
// isn't currently safe to go all the way down to 1 due to algorithmic
// assumptions.)
static const double kSwitchThresh = k2p800 * k2p50 * (1LL << 40);
static const double kLnSwitchThresh = 890.0 * kLn2;
static const double kJumpThresh = 314.0; // chosen to guarantee base_prob < kSwitchThresh in non-jumping case
// static const double kSwitchThresh = k2p100 * k2p50;
// static const double kLnSwitchThresh = 150.0 * kLn2;

// Since each Fisher's exact test contigency table has rows and columns, we use
// 'line' instead of 'row' to refer to the set of tables with 3rd column held
// constant.
BoolErr Fisher23LnFirstLine(int32_t obs_m11, int32_t obs_m12, int32_t obs_m21, int32_t obs_m22, double* tailp_ptr, dd_real* starting_lnprobv_ddr_ptr, int32_t* tie_ct_ptr, double* orig_base_probl_ptr, double* orig_base_lnprobl_ptr, double* orig_base_epsl_ptr, double* orig_base_probr_ptr, double* orig_base_lnprobr_ptr, double* orig_base_epsr_ptr, double* orig_saved_l11_ptr, double* orig_saved_l12_ptr, double* orig_saved_l21_ptr, double* orig_saved_l22_ptr, double* orig_saved_r11_ptr, double* orig_saved_r12_ptr, double* orig_saved_r21_ptr, double* orig_saved_r22_ptr) {
  // possible todo: have this and Fisher22LnP() call a shared function
  double m11 = obs_m11;
  double m12 = obs_m12;
  double m21 = obs_m21;
  double m22 = obs_m22;
  *starting_lnprobv_ddr_ptr =
    ddr_negate(ddr_add4_lfacts(m11, m12, m21, m22));
  double lastp = 1;
  int32_t tie_ct = 1;
  double tailp = 1;
  // No guaranteed normalization beyond m11+m21 >= m12+m22.  In particular,
  // m11 < m22 is possible.

  const double m1x = m11 + m12;
  const double m2x = m21 + m22;
  const double mx2 = m12 + m22;
  const double mxx = m1x + m2x;
  const double modal_m22 = m2x * mx2 / mxx;
  const double delta = m22 - modal_m22;
  // As with PLINK 1.9 fisher23(), "left tail" corresponds to small m22, and
  // "right tail" corresponds to large m22.
  if (delta > 0) {
    // Starting m22 is beginning of right tail.
    *orig_base_probr_ptr = 1;
    *orig_base_epsr_ptr = 0;
    *orig_saved_r11_ptr = m11;
    *orig_saved_r12_ptr = m12;
    *orig_saved_r21_ptr = m21;
    *orig_saved_r22_ptr = m22;
    // Iterate outward to floating-point precision limit.
    while (1) {
      m11 += 1;
      m22 += 1;
      lastp *= m12 * m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
    }
    if (delta < kJumpThresh) {
      // Jump back to starting table, and iterate inward until we find the
      // start of the other tail.
      lastp = 1;
      m11 = obs_m11;
      m12 = obs_m12;
      m21 = obs_m21;
      m22 = obs_m22;
      double one_plus_scaled_eps = 1;
      while (1) {
        const double m11_x_m22 = m11 * m22;
        // Don't want to wait until lastp becomes 0, since we want m11 >= 0 and
        // m22 >= 0 guaranteed when we save m11 and m22 on loop exit.
        if (m11_x_m22 == 0) {
          break;
        }
        m12 += 1;
        m21 += 1;
        lastp *= m11_x_m22 / (m12 * m21);
        m11 -= 1;
        m22 -= 1;
        one_plus_scaled_eps += 2 * k2m52;
        if (lastp < one_plus_scaled_eps) {
          if (lastp <= 2 - one_plus_scaled_eps) {
            tailp += lastp;
            break;
          }
          // Near-tie.  True value of lastp can be greater than, equal to, or
          // less than 1.
          const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
          intptr_t cmp_result;
          if (unlikely(Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, starting_lnprobv_ddr_ptr, &cmp_result, &lastp))) {
            return 1;
          }
          one_plus_scaled_eps = 1 + 3 * k2m52;
          if (cmp_result <= 0) {
            tailp += lastp;
            tie_ct += (cmp_result == 0);
            break;
          }
        }
      }
      *orig_saved_l11_ptr = m11;
      *orig_saved_l12_ptr = m12;
      *orig_saved_l21_ptr = m21;
      *orig_saved_l22_ptr = m22;
      *orig_base_probl_ptr = lastp;
      *orig_base_epsl_ptr = one_plus_scaled_eps - 1;
    } else {
      // Jump to other tail.  Round down and clamp.
      const double m22_minus_m11 = m22 - m11; // usually but not always nonpositive
      m22 = S_CAST(int32_t, modal_m22 - delta);
      const double min_m22 = MAXV(0, m22_minus_m11);
      if (m22 < min_m22) {
        m22 = min_m22;
      }
      while (1) {
        m11 = m22 - m22_minus_m11;
        m12 = mx2 - m22;
        m21 = m2x - m22;
        const dd_real lnprobv_ddr =
          ddr_negate(ddr_add4_lfacts(m11, m12, m21, m22));
        const dd_real lnprob_diff_ddr = ddr_sub(lnprobv_ddr, *starting_lnprobv_ddr_ptr);
        const double lnprob_diff = lnprob_diff_ddr.x[0];
        if (lnprob_diff >= k2m53) {
          if (m22 == min_m22) {
            // All tables on this tail have higher likelihood than the starting
            // table.  Exit.
            *tie_ct_ptr = 1;
            *orig_saved_l11_ptr = m11;
            *orig_saved_l12_ptr = m12;
            *orig_saved_l21_ptr = m21;
            *orig_saved_l22_ptr = m22;
            *tailp_ptr = tailp;
            if (lnprob_diff < kLnSwitchThresh) {
              *orig_base_probl_ptr = ddr_exp(lnprob_diff_ddr).x[0];
              *orig_base_epsl_ptr = k2m52;
            } else {
              *orig_base_probl_ptr = 0;
              *orig_base_lnprobl_ptr = lnprob_diff;
              *orig_base_epsl_ptr = ceil_smalleps(lnprob_diff) * k2m52;
            }
            return 0;
          }
          const double ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
          // Round up, to guarantee that we make progress.
          // This may overshoot.  But the function is guaraneed to terminate
          // because we never overshoot (and we do always make progress on each
          // step) once we're on the other side.
          m22 -= ceil_smalleps(lnprob_diff / ll_deriv);
          if (m22 < min_m22) {
            m22 = min_m22;
          }
        } else if (lnprob_diff > -62 * kLn2) {
          lastp = exp(lnprob_diff);
          break;
        } else {
          const double ll_deriv = log(m12 * m21 / ((m11 + 1) * (m22 + 1)));
          // Round down, to guarantee we don't overshoot.
          // (lnprob_diff is negative and ll_deriv is positive.)
          // We're guaranteed to make progress, since lnprob_diff <=
          // -62 * log(2) and m12 * m21 >= 1.
          m22 -= S_CAST(int64_t, lnprob_diff / ll_deriv);
        }
      }
      // Sum toward center, until lastp >= 1.
      double one_minus_scaled_eps = 1 - 3 * k2m52;
      const double lastp_tail = lastp;
      const double m11_tail = m11;
      const double m12_tail = m12;
      const double m21_tail = m21;
      const double m22_tail = m22;
      while (lastp <= one_minus_scaled_eps) {
        tailp += lastp;
        m11 += 1;
        m22 += 1;
        lastp *= m12 * m21 / (m11 * m22);
        m12 -= 1;
        m21 -= 1;
        one_minus_scaled_eps -= 2 * k2m52;
      }
      if (lastp < 2 - one_minus_scaled_eps) {
        const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
        intptr_t cmp_result;
        if (unlikely(Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, starting_lnprobv_ddr_ptr, &cmp_result, &lastp))) {
          return 1;
        }
        one_minus_scaled_eps = 1 - 3 * k2m52;
        if (cmp_result <= 0) {
          tailp += lastp;
          tie_ct += (cmp_result == 0);
        }
      }
      *orig_saved_l11_ptr = m11;
      *orig_saved_l12_ptr = m12;
      *orig_saved_l21_ptr = m21;
      *orig_saved_l22_ptr = m22;
      *orig_base_probl_ptr = lastp;
      *orig_base_epsl_ptr = 1 - one_minus_scaled_eps;
      lastp = lastp_tail;
      m11 = m11_tail;
      m12 = m12_tail;
      m21 = m21_tail;
      m22 = m22_tail;
    }
    *tie_ct_ptr = tie_ct;
    // Sum away from center, until sums stop changing.
    while (1) {
      m12 += 1;
      m21 += 1;
      lastp *= m11 * m22 / (m12 * m21);
      const double preaddp = tailp;
      tailp += lastp;
      if (tailp == preaddp) {
        break;
      }
      m11 -= 1;
      m22 -= 1;
    }
    *tailp_ptr = tailp;
    return 0;
  }
  *orig_base_probl_ptr = 1;
  *orig_base_epsl_ptr = 0;
  *orig_saved_l11_ptr = obs_m11;
  *orig_saved_l12_ptr = obs_m12;
  *orig_saved_l21_ptr = obs_m21;
  *orig_saved_l22_ptr = obs_m22;
  while (1) {
    m12 += 1;
    m21 += 1;
    lastp *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
  }
  if (delta > -kJumpThresh) {
    lastp = 1;
    m11 = obs_m11;
    m12 = obs_m12;
    m21 = obs_m21;
    m22 = obs_m22;
    double one_plus_scaled_eps = 1;
    while (1) {
      const double m12_x_m21 = m12 * m21;
      if (m12_x_m21 == 0) {
        break;
      }
      m11 += 1;
      m22 += 1;
      lastp *= m12_x_m21 / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      one_plus_scaled_eps += 2 * k2m52;
      if (lastp < one_plus_scaled_eps) {
        if (lastp <= 2 - one_plus_scaled_eps) {
          tailp += lastp;
          break;
        }
        const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
        intptr_t cmp_result;
        if (unlikely(Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, starting_lnprobv_ddr_ptr, &cmp_result, &lastp))) {
          return 1;
        }
        one_plus_scaled_eps = 1 + 3 * k2m52;
        if (cmp_result <= 0) {
          tailp += lastp;
          tie_ct += (cmp_result == 0);
          break;
        }
      }
    }
    *orig_saved_r11_ptr = m11;
    *orig_saved_r12_ptr = m12;
    *orig_saved_r21_ptr = m21;
    *orig_saved_r22_ptr = m22;
    *orig_base_probr_ptr = lastp;
    *orig_base_epsr_ptr = one_plus_scaled_eps - 1;
  } else {
    // Jump to other tail.
    const double m21_minus_m12 = m21 - m12; // usually but not always nonpositive
    const double mx1 = obs_m11 + obs_m21;
    const double modal_m21 = m2x - modal_m22;
    m21 = S_CAST(int32_t, modal_m21 + delta);
    const double min_m21 = MAXV(0, m21_minus_m12);
    if (m21 < min_m21) {
      m21 = min_m21;
    }
    while (1) {
      m11 = mx1 - m21;
      m12 = m21 - m21_minus_m12;
      m22 = m2x - m21;
      const dd_real lnprobv_ddr =
        ddr_negate(ddr_add4_lfacts(m11, m12, m21, m22));
      const dd_real lnprob_diff_ddr = ddr_sub(lnprobv_ddr, *starting_lnprobv_ddr_ptr);
      const double lnprob_diff = lnprob_diff_ddr.x[0];
      if (lnprob_diff >= k2m53) {
        if (m21 == min_m21) {
          // All tables on this tail have higher likelihood than the starting
          // table.  Exit.
          *tie_ct_ptr = 1;
          *orig_saved_r11_ptr = m11;
          *orig_saved_r12_ptr = m12;
          *orig_saved_r21_ptr = m21;
          *orig_saved_r22_ptr = m22;
          *tailp_ptr = tailp;
          if (lnprob_diff < kLnSwitchThresh) {
            *orig_base_probr_ptr = ddr_exp(lnprob_diff_ddr).x[0];
            *orig_base_epsr_ptr = k2m52;
          } else {
            *orig_base_probr_ptr = 0;
            *orig_base_lnprobr_ptr = lnprob_diff;
            *orig_base_epsr_ptr = ceil_smalleps(lnprob_diff) * k2m52;
          }
          return 0;
        }
        // Derivative is w.r.t. m21.
        const double ll_deriv = log((m11 + 1) * (m22 + 1) / (m12 * m21));
        m21 -= ceil_smalleps(lnprob_diff / ll_deriv);
        if (m21 < min_m21) {
          m21 = min_m21;
        }
      } else if (lnprob_diff > -62 * kLn2) {
        lastp = exp(lnprob_diff);
        break;
      } else {
        const double ll_deriv = log(m11 * m22 / ((m12 + 1) * (m21 + 1)));
        m21 += S_CAST(int64_t, -lnprob_diff / ll_deriv);
      }
    }
    // Sum toward center, until lastp >= 1.
    double one_minus_scaled_eps = 1 - 3 * k2m52;
    const double lastp_tail = lastp;
    const double m11_tail = m11;
    const double m12_tail = m12;
    const double m21_tail = m21;
    const double m22_tail = m22;
    while (lastp <= one_minus_scaled_eps) {
      tailp += lastp;
      m12 += 1;
      m21 += 1;
      lastp *= m11 * m22 / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      one_minus_scaled_eps -= 2 * k2m52;
    }
    if (lastp < 2 - one_minus_scaled_eps) {
      const int32_t m22_incr = S_CAST(int32_t, m22) - obs_m22;
      intptr_t cmp_result;
      if (unlikely(Fisher22Compare(obs_m11, obs_m12, obs_m21, obs_m22, m22_incr, starting_lnprobv_ddr_ptr, &cmp_result, &lastp))) {
        return 1;
      }
      one_minus_scaled_eps = 1 - 3 * k2m52;
      if (cmp_result <= 0) {
        tailp += lastp;
        tie_ct += (cmp_result == 0);
      }
    }
    *orig_saved_r11_ptr = m11;
    *orig_saved_r12_ptr = m12;
    *orig_saved_r21_ptr = m21;
    *orig_saved_r22_ptr = m22;
    *orig_base_probr_ptr = lastp;
    *orig_base_epsr_ptr = 1 - one_minus_scaled_eps;
    lastp = lastp_tail;
    m11 = m11_tail;
    m12 = m12_tail;
    m21 = m21_tail;
    m22 = m22_tail;
  }
  *tie_ct_ptr = tie_ct;
  // Sum away from center, until sums stop changing.
  while (1) {
    m11 += 1;
    m22 += 1;
    lastp *= m12 * m21 / (m11 * m22);
    const double preaddp = tailp;
    tailp += lastp;
    if (tailp == preaddp) {
      break;
    }
    m12 -= 1;
    m21 -= 1;
  }
  *tailp_ptr = tailp;
  return 0;
}

BoolErr Fisher23Compare(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m13, uint32_t obs_m21, uint32_t obs_m22, uint32_t obs_m23, uint32_t cur_m11, uint32_t cur_m12, dd_real* neg_numer_ddr_ptr, intptr_t* cmp_resultp, double* dbl_ptr) {
  // Fisher 2x3 likelihood is
  //
  //   (m11+m12+m13)! (m21+m22+m23)! (m11+m21)! (m12+m22)! (m13+m23)!
  //   --------------------------------------------------------------
  //      m11! m12! m13! m21! m22! m23! (m11+m12+m13+m21+m22+m23)!
  //
  // so likelihood ratio of interest is
  //
  //   obs_m11! obs_m12! obs_m13! obs_m21! obs_m22! obs_m23!
  //   -----------------------------------------------------
  //   cur_m11! cur_m12! cur_m13! cur_m21! cur_m22! cur_m23!
  //
  // where we infer the remaining cur_ values from the fixed row/column sums.
  const uint32_t cur_m13 = obs_m11 + obs_m12 + obs_m13 - cur_m11 - cur_m12;
  const uint32_t cur_m21 = obs_m11 + obs_m21 - cur_m11;
  const uint32_t cur_m22 = obs_m12 + obs_m22 - cur_m12;
  const uint32_t cur_m23 = obs_m13 + obs_m23 - cur_m13;
  uint32_t numer_factorial_args[6];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m13;
  numer_factorial_args[3] = obs_m21;
  numer_factorial_args[4] = obs_m22;
  numer_factorial_args[5] = obs_m23;
  uint32_t denom_factorial_args[6];
  denom_factorial_args[0] = cur_m11;
  denom_factorial_args[1] = cur_m12;
  denom_factorial_args[2] = cur_m13;
  denom_factorial_args[3] = cur_m21;
  denom_factorial_args[4] = cur_m22;
  denom_factorial_args[5] = cur_m23;

  mp_limb_t* gmp_wkspace = nullptr;
  uintptr_t gmp_wkspace_limb_ct = 0;
  BoolErr reterr = CompareFactorialProducts(6, 0, 0, numer_factorial_args, denom_factorial_args, neg_numer_ddr_ptr, &gmp_wkspace, &gmp_wkspace_limb_ct, cmp_resultp, dbl_ptr);
  free_cond(gmp_wkspace);
  return reterr;
}

// Internally, this assumes we are iterating down the small-m22 tail.  Pass in
// saved_m1lo = m12, saved_m1hi = m11, saved_m2hi = m22, saved_m2lo = m21 to
// iterate down the large-m22 tail; this also requires swapping the
// corresponding obs_ values.
BoolErr Fisher23LnPTailsum(dd_real starting_lnprobv_ddr, uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m13, uint32_t obs_m21, uint32_t obs_m22, uint32_t obs_m23, double* base_probp, double* base_lnprobp, double* base_epsp, double* saved_m1lop, double* saved_m1hip, double* saved_m2hip, double* saved_m2lop, int32_t* tie_ctp, double* totalp, uint32_t* center_is_emptyp) {
  double total = 0;
  double lastp = *base_probp;
  double cur_eps = *base_epsp;
  double m11 = *saved_m1lop;
  double m12 = *saved_m1hip;
  double m21 = *saved_m2hip;
  double m22 = *saved_m2lop;
  // identify beginning (center-facing side) of tail
  if (lastp == 0.0) {
    double last_lnp = *base_lnprobp;
    if (last_lnp >= kLnSwitchThresh) {
      while (1) {
        const double prev_numer = m11 * m22;
        if (prev_numer == 0) {
          // lowest-likelihood table on this side is still too probable
          *base_lnprobp = last_lnp;
          *base_epsp = cur_eps;
          *center_is_emptyp = 0;
          *saved_m1lop = m11;
          *saved_m1hip = m12;
          *saved_m2hip = m21;
          *saved_m2lop = m22;
          return 0;
        }
        m12 += 1;
        m21 += 1;
        const double lnprob_incr = log(prev_numer / (m12 * m21));
        last_lnp += lnprob_incr;
        m11 -= 1;
        m21 -= 1;
        cur_eps += 3 * k2m52;
        if (lnprob_incr <= -2) {
          cur_eps += (trunc(-lnprob_incr) - 1) * k2m52;
        }
        // no risk of last_lnp dropping below cur_eps
      }
    }
    lastp = exp(last_lnp);
  }
  double m11_tail;
  double m12_tail;
  double m21_tail;
  double m22_tail;
  if (lastp >= 1 + cur_eps) {
    while (1) {
      const double prev_numer = m11 * m22;
      if (prev_numer == 0) {
        // lowest-likelihood table on this side is still too probable
        *center_is_emptyp = 0;
        *saved_m1lop = m11;
        *saved_m1hip = m12;
        *saved_m2hip = m21;
        *saved_m2lop = m22;
        if (lastp < kSwitchThresh) {
          *base_probp = lastp;
          *base_epsp = cur_eps;
        } else {
          *base_probp = 0;
          const double last_lnp = log(lastp);
          *base_lnprobp = last_lnp;
          *base_epsp = cur_eps + ceil_smalleps(last_lnp) * k2m52;
        }
        return 0;
      }
      m12 += 1;
      m21 += 1;
      lastp *= prev_numer / (m12 * m21);
      m11 -= 1;
      m22 -= 1;
      cur_eps += 2 * k2m52;
      if (lastp < 1 + cur_eps) {
        if (lastp <= 1 - cur_eps) {
          break;
        }
        intptr_t cmp_result;
        if (unlikely(Fisher23Compare(obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, S_CAST(int32_t, m11), S_CAST(int32_t, m12), &starting_lnprobv_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        cur_eps = 3 * k2m52;
        if (cmp_result <= 0) {
          *tie_ctp += (cmp_result == 0);
          break;
        }
      }
    }
    *base_probp = lastp;
    total = lastp;
    m11_tail = m11;
    m12_tail = m12;
    m21_tail = m21;
    m22_tail = m22;
  } else {
    m11_tail = m11;
    m12_tail = m12;
    m21_tail = m21;
    m22_tail = m22;
    const double lastp_tail = lastp;
    uint32_t tie_ct_incr = 0;
    while (1) {
      if (lastp > 1 - cur_eps) {
        if (lastp >= 1 + cur_eps) {
          break;
        }
        intptr_t cmp_result;
        if (unlikely(Fisher23Compare(obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, S_CAST(int32_t, m11), S_CAST(int32_t, m12), &starting_lnprobv_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        cur_eps = 3 * k2m52;
        if (cmp_result > 0) {
          break;
        }
        tie_ct_incr += (cmp_result == 0);
      }
      total += lastp;
      m11 += 1;
      m22 += 1;
      const double prob_mult = m12 * m21 / (m11 * m22);
      if (prob_mult <= 1 + 2 * k2m52) {
        const int64_t m11_i = S_CAST(int32_t, m11);
        const int64_t m12_i = S_CAST(int32_t, m12);
        const int64_t m21_i = S_CAST(int32_t, m21);
        const int64_t m22_i = S_CAST(int32_t, m22);
        const int64_t prob_mult_numer = m12_i * m21_i;
        const int64_t prob_mult_denom = m11_i * m22_i;
        if (prob_mult_numer <= prob_mult_denom) {
          if (tie_ct_incr) {
            *tie_ctp += tie_ct_incr + (prob_mult_numer == prob_mult_denom);
          }
          *center_is_emptyp = 1;
          return 0;
        }
      }
      lastp *= prob_mult;
      m12 -= 1;
      m21 -= 1;
      cur_eps += 2 * k2m52;
    }
    *base_probp = lastp;
    *tie_ctp += tie_ct_incr;
    lastp = lastp_tail;
  }
  *base_epsp = cur_eps;
  *saved_m1lop = m11;
  *saved_m1hip = m12;
  *saved_m2hip = m21;
  *saved_m2lop = m22;
  // {m11,m12,m21,m22}_tail is now at a table with relative-likelihood <= 1,
  // and 'total' is the sum of all relative-likelihoods for tables at least as
  // close to the center.

  // sum tail to floating-point precision limit
  m11 = m11_tail;
  m12 = m12_tail;
  m21 = m21_tail;
  m22 = m22_tail;
  while (1) {
    m12 += 1;
    m21 += 1;
    lastp *= m11 * m22 / (m12 * m21);
    m11 -= 1;
    m22 -= 1;
    const double preaddp = total;
    total += lastp;
    if (total == preaddp) {
      break;
    }
  }
  *totalp = total;
  *center_is_emptyp = 0;
  return 0;
}

// obs_m11 + obs_m12 + obs_m13 + obs_m21 + obs_m22 + obs_m23 assumed to be
// <2^31.
BoolErr Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp, double* resultp) {
  // Normalize: m11 + m21 >= m12 + m22 >= m13 + m23,
  //            m13 * (m21 + m22) <= m23 * (m11 + m12)
  // Note that columns are reversed from PLINK 1.9 fisher23().
  {
    int32_t mx1i = obs_m11 + obs_m21;
    {
      int32_t mx2i = obs_m12 + obs_m22;
      if (mx1i < mx2i) {
        swap_i32(&obs_m11, &obs_m12);
        swap_i32(&obs_m21, &obs_m22);
        swap_i32(&mx1i, &mx2i);
      }
      {
        int32_t mx3i = obs_m13 + obs_m23;
        if (mx2i < mx3i) {
          swap_i32(&obs_m12, &obs_m13);
          swap_i32(&obs_m22, &obs_m23);
          swap_i32(&mx2i, &mx3i);
        }
        if (mx3i == 0) {
          return Fisher22LnP(obs_m11, obs_m12, obs_m21, obs_m22, midp, resultp);
        }
      }
      if (mx1i < mx2i) {
        swap_i32(&obs_m11, &obs_m12);
        swap_i32(&obs_m21, &obs_m22);
        mx1i = mx2i;
      }
    }
    if (S_CAST(int64_t, obs_m13) * (obs_m21 + obs_m22) > S_CAST(int64_t, obs_m23) * (obs_m11 + obs_m12)) {
      swap_i32(&obs_m11, &obs_m21);
      swap_i32(&obs_m12, &obs_m22);
      swap_i32(&obs_m13, &obs_m23);
    }
  }
  double orig_base_lnprobl = 0;
  double orig_base_lnprobr = 0;
  double orig_base_probl;
  double orig_base_epsl;
  double orig_base_probr;
  double orig_base_epsr;
  double tailp;
  dd_real starting_lnprobv_ddr;
  int32_t tie_ct;
  double orig_saved_l11;
  double orig_saved_l12;
  double orig_saved_l21;
  double orig_saved_l22;
  double orig_saved_r11;
  double orig_saved_r12;
  double orig_saved_r21;
  double orig_saved_r22;
  if (unlikely(Fisher23LnFirstLine(obs_m11, obs_m12, obs_m21, obs_m22, &tailp, &starting_lnprobv_ddr, &tie_ct, &orig_base_probl, &orig_base_lnprobl, &orig_base_epsl, &orig_base_probr, &orig_base_lnprobr, &orig_base_epsr, &orig_saved_l11, &orig_saved_l12, &orig_saved_l21, &orig_saved_l22, &orig_saved_r11, &orig_saved_r12, &orig_saved_r21, &orig_saved_r22))) {
    return 1;
  }

  // Returned starting_lnprobv_ddr is for a single line, corresponding to
  // 1 / (obs_m11! obs_m12! obs_m21! obs_m22!).
  // Include m13! and m23! before we iterate over other lines.
  starting_lnprobv_ddr =
    ddr_sub(starting_lnprobv_ddr,
            ddr_add_lfacts(obs_m13, obs_m23));

  // Other log-factorial expressions we want:
  //
  // * lnprobf, so we can add it to starting_lnprobv at the end and convert
  //   tailp (which is a multiple of starting_prob) into the final
  //   log-[mid]p-value.
  //
  //     (m11+m12+m13)! (m21+m22+m23)! (m11+m21)! (m12+m22)! (m13+m23)!
  //     --------------------------------------------------------------
  //                       (m11+m12+m13+m21+m22+m23)!
  //
  // * Likelihood sum over a 3rd-column-constant line is
  //
  //     (m11+m12+m13)! (m21+m22+m23)! (m13+m23)!    (m11+m12+m21+m22)!
  //     ---------------------------------------- * ---------------------
  //       (m11+m12+m13+m21+m22+m23)! m13! m23!     (m11+m12)! (m21+m22)!
  //
  //   i.e. the 2x2 likelihood after merging columns 1 and 2.
  //
  //   Converting to starting_prob units:
  //
  //                    1                           (m11+m12+m21+m22)!
  //     ------------------------------- * ------------------------------------
  //     m13! m23! (m11+m12)! (m21+m22)!   (m11+m21)! (m12+m22)! starting_probv
  const double m1x = obs_m11 + obs_m12 + obs_m13;
  const double m2x = obs_m21 + obs_m22 + obs_m23;
  const double mx1 = obs_m11 + obs_m21;
  const double mx2 = obs_m12 + obs_m22;
  const double mx3 = obs_m13 + obs_m23;
  const double mxx = m1x + m2x;
  const dd_real lnprobf_ddr =
    ddr_sub(ddr_add5_lfacts(m1x, m2x, mx1, mx2, mx3),
            ddr_lfact(mxx));
  const dd_real line_relative_lnprobf_ddr =
    ddr_sub(ddr_lfact(mx1 + mx2),
            ddr_add(ddr_add_lfacts(mx1, mx2), starting_lnprobv_ddr));

  // tailp, orig_base_probl, and orig_base_probr are in starting_prob units.
  for (uint32_t m13_decreasing = 0; m13_decreasing != 2; ++m13_decreasing) {
    double m13 = obs_m13;
    double m23 = obs_m23;
    // left = small m11 and m22; right = large m11 and m22.
    double l11 = orig_saved_l11;
    double l12 = orig_saved_l12;
    double l21 = orig_saved_l21;
    double l22 = orig_saved_l22;
    double r11 = orig_saved_r11;
    double r12 = orig_saved_r12;
    double r21 = orig_saved_r21;
    double r22 = orig_saved_r22;
    double base_probl = orig_base_probl;
    double base_lnprobl = orig_base_lnprobl;
    double base_epsl = orig_base_epsl;
    double base_probr = orig_base_probr;
    double base_lnprobr = orig_base_lnprobr;
    double base_epsr = orig_base_epsr;
    int32_t m13_decr_last;
    if (m13_decreasing) {
      m13_decr_last = MINV(obs_m13, obs_m21 + obs_m22);
    } else {
      m13_decr_last = -MINV(obs_m23, obs_m11 + obs_m12);
    }
    for (int32_t m13_decr = 0; m13_decr != m13_decr_last; ) {
      m13_decr += m13_decreasing * 2 - 1;
      double prob_mult;
      if (m13_decreasing) {
        m23 += 1;
        // Need to be careful to not move l{11,12,21,22} to the right of the
        // mode, or r{11,12,21,22} to the left of it.
        if (l22) {
          l12 += 1;
          prob_mult = m13 * l22 / (m23 * l12);
          l22 -= 1;
        } else {
          l11 += 1;
          prob_mult = m13 * l21 / (m23 * l11);
          l21 -= 1;
        }
        m13 -= 1;
      } else {
        m13 += 1;
        if (l11) {
          l21 += 1;
          prob_mult = m23 * l11 / (m13 * l21);
          l11 -= 1;
        } else {
          l22 += 1;
          prob_mult = m23 * l12 / (m13 * l22);
          l12 -= 1;
        }
        m23 -= 1;
      }
      if (base_probl == 0.0) {
        const double lnprob_incr = log(prob_mult);
        base_lnprobl += lnprob_incr;
        base_epsl += 3 * k2m52;
        if (fabs(lnprob_incr) >= 2) {
          base_epsl += (trunc(fabs(lnprob_incr)) - 1) * k2m52;
        }
      } else {
        base_probl *= prob_mult;
        base_epsl += 2 * k2m52;
      }
      double tail_incr1 = 0.0;
      uint32_t center_is_empty;
      if (unlikely(Fisher23LnPTailsum(starting_lnprobv_ddr, obs_m11, obs_m12, obs_m13, obs_m21, obs_m22, obs_m23, &base_probl, &base_lnprobl, &base_epsl, &l11, &l12, &l21, &l22, &tie_ct, &tail_incr1, &center_is_empty))) {
        return 1;
      }
      if (center_is_empty) {
        // All tables in this line, and all subsequent lines, are less probable
        // than the initial table.
        double m11_12 = m1x - m13;
        double m21_22 = m2x - m23;
        const dd_real line_adj_ddr = ddr_add4_lfacts(m13, m23, m11_12, m21_22);
        double line_prob = exp(ddr_sub(line_relative_lnprobf_ddr, line_adj_ddr).x[0]);
        if (m13_decreasing) {
          while (1) {
            const double preaddp = tailp;
            tailp += line_prob;
            if (tailp == preaddp) {
              break;
            }
            m11_12 += 1;
            m23 += 1;
            line_prob *= m13 * m21_22 / (m23 * m11_12);
            m13 -= 1;
            m21_22 -= 1;
          }
        } else {
          while (1) {
            const double preaddp = tailp;
            tailp += line_prob;
            if (tailp == preaddp) {
              break;
            }
            m13 += 1;
            m21_22 += 1;
            line_prob *= m11_12 * m23 / (m13 * m21_22);
            m11_12 -= 1;
            m23 -= 1;
          }
        }
        break;
      }
      tailp += tail_incr1;
      if (m13_decreasing) {
        const double prev_m13 = m13 + 1;
        if (r21) {
          r11 += 1;
          prob_mult = prev_m13 * r21 / (m23 * r11);
          r21 -= 1;
        } else {
          r12 += 1;
          prob_mult = prev_m13 * r22 / (m23 * r12);
          r22 -= 1;
        }
      } else {
        const double prev_m23 = m23 + 1;
        if (r12) {
          r22 += 1;
          prob_mult = prev_m23 * r12 / (m13 * r22);
          r12 -= 1;
        } else {
          r21 += 1;
          prob_mult = prev_m23 * r11 / (m13 * r21);
          r11 -= 1;
        }
      }
      if (base_probr == 0.0) {
        const double lnprob_incr = log(prob_mult);
        base_lnprobr += lnprob_incr;
        base_epsr += 3 * k2m52;
        if (fabs(lnprob_incr) >= 2) {
          base_epsr += (trunc(fabs(lnprob_incr)) - 1) * k2m52;
        }
      } else {
        base_probr *= prob_mult;
        base_epsr += 2 * k2m52;
      }
      double tail_incr2 = 0.0;
      if (unlikely(Fisher23LnPTailsum(starting_lnprobv_ddr, obs_m12, obs_m11, obs_m13, obs_m22, obs_m21, obs_m23, &base_probr, &base_lnprobr, &base_epsr, &r12, &r11, &r22, &r21, &tie_ct, &tail_incr2, &center_is_empty))) {
        return 1;
      }
      tailp += tail_incr2;
    }
  }
  const double starting_lnprob =
    ddr_add(lnprobf_ddr, starting_lnprobv_ddr).x[0];
  if (midp) {
    tailp -= tie_ct * 0.5;
  }
  const double result = log(tailp) + starting_lnprob;
  *resultp = result;
  if (result > -k2m35) {
    // true p-value should always be 1 here
    // (possible todo: check boundary cases with total near 2^31)
    *resultp = 0;
  }
  return 0;
}

#ifdef __cplusplus
}
#endif
