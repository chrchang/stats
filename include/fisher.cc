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

// *cmp_resultp is set to positive value if m11 = obs_m11 + m11_incr has higher
// likelihood than m11 = obs_m11, 0 if identical likelihood, and negative value
// if lower likelihood.
// Error is returned iff malloc fails.
// If neg_numer_ddr has not been computed yet, set its
// x[0] to DBL_MAX; it will be filled in if necessary.
BoolErr FisherCompare(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m21, uint32_t obs_m22, int32_t m11_incr, dd_real* neg_numer_ddr_ptr, intptr_t* cmp_resultp, double* dbl_ptr) {
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
  // where j=m11_incr.
  uint32_t numer_factorial_args[4];
  numer_factorial_args[0] = obs_m11;
  numer_factorial_args[1] = obs_m12;
  numer_factorial_args[2] = obs_m21;
  numer_factorial_args[3] = obs_m22;
  uint32_t denom_factorial_args[4];
  denom_factorial_args[0] = obs_m11 + m11_incr;
  denom_factorial_args[0] = obs_m12 - m11_incr;
  denom_factorial_args[0] = obs_m21 - m11_incr;
  denom_factorial_args[0] = obs_m22 + m11_incr;

  mp_limb_t* gmp_wkspace = nullptr;
  uintptr_t gmp_wkspace_limb_ct = 0;
  BoolErr reterr = CompareFactorialProducts(4, 0, 0, numer_factorial_args, denom_factorial_args, neg_numer_ddr_ptr, &gmp_wkspace, &gmp_wkspace_limb_ct, cmp_resultp, dbl_ptr);
  free_cond(gmp_wkspace);
  return reterr;
}

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^31.
BoolErr Fisher22LnP(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m21, uint32_t obs_m22, uint32_t midp, double* resultp) {
  // Normalize.
  if (obs_m11 > obs_m22) {
    swap_u32(&obs_m11, &obs_m22);
  }
  if (obs_m12 > obs_m21) {
    swap_u32(&obs_m12, &obs_m21);
  }
  if (S_CAST(uint64_t, obs_m11) * obs_m22 > S_CAST(uint64_t, obs_m12) * obs_m21) {
    swap_u32(&obs_m11, &obs_m12);
    swap_u32(&obs_m21, &obs_m22);
  }
  if (!midp) {
    // Fast path for p=1.
    if (S_CAST(uint64_t, obs_m11 + 1) * (obs_m22 + 1) == S_CAST(uint64_t, obs_m12) * obs_m21) {
      *resultp = 0;
      return 0;
    }
  }
  // Iterate outward to floating-point precision limit.
  const double obs_m11d = u31tod(obs_m11);
  const double obs_m12d = u31tod(obs_m12);
  const double obs_m21d = u31tod(obs_m21);
  const double obs_m22d = u31tod(obs_m22);
  double m11 = obs_m11d;
  double m12 = obs_m12d;
  double m21 = obs_m21d;
  double m22 = obs_m22d;
  double lastp = 1;
  double tailp = 1;
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
  int32_t tie_ct = 1;
  if ((S_CAST(int64_t, obs_m11) + 172) * (S_CAST(int64_t, obs_m22) + 172) >=
      (S_CAST(int64_t, obs_m12) - 172) * (S_CAST(int64_t, obs_m21) - 172)) {
    lastp = 1;
    m11 = obs_m11d;
    m12 = obs_m12d;
    m21 = obs_m21d;
    m22 = obs_m22d;
    dd_real starting_lnprob_other_component_ddr = {{DBL_MAX, 0.0}};
    double centerp = 0;
    while (1) {
      m11 += 1;
      m22 += 1;
      lastp *= (m12 * m21) / (m11 * m22);
      m12 -= 1;
      m21 -= 1;
      // Number of center contingency tables is maximized with obs_m11 = 0,
      // modal_m11 = 172, other values large.
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
        const intptr_t m11_incr = S_CAST(intptr_t, m11) - obs_m11;
        intptr_t cmp_result;
        if (unlikely(FisherCompare(obs_m11, obs_m12, obs_m21, obs_m22, m11_incr, &starting_lnprob_other_component_ddr, &cmp_result, &lastp))) {
          return 1;
        }
        if (cmp_result <= 0) {
          tailp += lastp;
          tie_ct += (cmp_result == 0);
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
    const double denom = tailp + centerp;
    if (midp) {
      tailp -= S_CAST(double, tie_ct) * 0.5;
    }
    *resultp = log(tailp / denom);
    return 0;
  }
  dd_real starting_lnprob_other_component_ddr =
    ddr_negate(ddr_add4_lfacts(obs_m11d, obs_m12d, obs_m21d, obs_m22d));

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
  // The current heuristic starts by reflecting (obs_m12 + m12) * 0.5 across
  // the mode, performing a full log-likelihood check at the nearest valid
  // point.  (It is convenient to focus on m12 here, since m12=0 corresponds to
  // the outermost table on this tail.)  Hopefully we find that we're in
  // (starting_lnprob - 62 * kLn2, starting_lnprob], so we're at or near a
  // table that actually contributes to the tail-sum.  (This window is chosen
  // to be wide enough to guarantee that at least one point falls inside when
  // obs_m11 + obs_m12 + obs_m21 + obs_m22 < 2^31.)
  //
  // If not, we jump again, using Newton's method.
  // If m12 is too high (i.e. current log-likelihood is too high), decreasing
  // m12 by 1 would multiply the likelihood by
  //   (m11 + 1) * (m22 + 1) / (m12 * m21)
  // If m12 is too low, increasing m12 by 1 would multiply the likelihood by
  //   (m12 + 1) * (m21 + 1) / (m11 * m22)
  // We use the negative-log of the first expression as the Newton's method
  // f'(x) when we're jumping to lower m12, and the log of the second
  // expression when we're jumping to higher m12.
  // f''(x) is always negative, so we can aim for starting_lnprob instead of
  // the middle of the interval.

  const double m21_minus_m12 = obs_m21d - obs_m12d;
  const double m1x = obs_m11d + obs_m12d;
  const double m2x = obs_m21d + obs_m22d;
  const double mx2 = obs_m12d + obs_m22d;
  const double mxx = m1x + m2x;
  {
    // x=modal_m12 satisfies
    //    (m1x - x) * (mx2 - x) = x * (x + m21_minus_m12)
    // -> (x - m1x) * (x - mx2) = x * (x + m21_minus_m12)
    // -> x^2 + x*(-m1x - mx2) + m1x*mx2 = x^2 + x*m21_minus_m12
    // -> m1x*mx2 = x*(m21_minus_m12 + m1x + mx2)
    // -> x = m1x*mx2 / mxx
    const double modal_m12 = m1x * mx2 / mxx;
    m12 = 2 * modal_m12 - (m12 + obs_m12d) * 0.5;
    // Round down (to guarantee we've actually moved to the other side of the
    // mode) and clamp.
    m12 = S_CAST(double, S_CAST(int32_t, m12));
    if (m12 < 0) {
      m12 = 0;
    }
  }
  const dd_real common_lnprob_component_ddr =
    ddr_sub(ddr_add4_lfacts(m1x, m2x, obs_m11d + obs_m21d, mx2),
            ddr_lfact(mxx));
  const double starting_lnprob = ddr_add(common_lnprob_component_ddr, starting_lnprob_other_component_ddr).x[0];
  while (1) {
    m11 = m1x - m12;
    m21 = m21_minus_m12 + m12;
    m22 = mx2 - m12;
    const dd_real lnprob_other_component_ddr =
      ddr_negate(ddr_add4_lfacts(m11, m12, m21, m22));
    const double lnprob_diff = ddr_sub(lnprob_other_component_ddr, starting_lnprob_other_component_ddr).x[0];
    // Could tighten this threshold further; I haven't performed a careful
    // error analysis yet but CompareFactorialProducts() includes a plausible
    // assumption that 2^{-60} is safe.  But code is correct as long as we're
    // guaranteed to enter the "lastp < 2 - one_minus_scaled_eps" branch for
    // positive lnprob_diff.
    if (lnprob_diff >= k2m53) {
      if (m12 == 0) {
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
      m12 -= 1 - S_CAST(int64_t, (1 - kSmallEpsilon) * lnprob_diff / ll_deriv);
      if (m12 < 0) {
        m12 = 0;
      }
    } else if (lnprob_diff > -62 * kLn2) {
      lastp = exp(lnprob_diff);
      break;
    } else {
      const double ll_deriv = log((m12 + 1) * (m21 + 1) / (m11 * m22));
      // Round down, to guarantee we don't overshoot.
      // We're guaranteed to make progress, since lnprob_diff <= -62 * log(2),
      // m11 * m22 < 2^62, and (m12 + 1) * (m21 + 1) >= 1.
      m12 += S_CAST(int64_t, lnprob_diff / ll_deriv);
    }
  }
  // Sum toward center, until lastp >= 1.
  // lastp should be accurate to 3 ULP as we enter this loop (max 1.5 ULP
  // observed error from exp, tiny bit over 0.5 from lnprob_diff), so near-tie
  // detection can use a tight epsilon here.
  double one_minus_scaled_eps = 1 - 3 * k2m52;
  double lastp_tail = lastp;
  double m11_center = m11;
  double m12_center = m12;
  double m21_center = m21;
  double m22_center = m22;
  while (lastp <= one_minus_scaled_eps) {
    tailp += lastp;
    m12_center += 1;
    m21_center += 1;
    lastp *= m11_center * m22_center / (m12_center * m21_center);
    m11_center -= 1;
    m22_center -= 1;
    one_minus_scaled_eps -= 2 * k2m52;
  }
  if (lastp < 2 - one_minus_scaled_eps) {
    const intptr_t m11_incr = S_CAST(intptr_t, m11) - obs_m11;
    intptr_t cmp_result;
    if (unlikely(FisherCompare(obs_m11, obs_m12, obs_m21, obs_m22, m11_incr, &starting_lnprob_other_component_ddr, &cmp_result, &lastp))) {
      return 1;
    }
    if (cmp_result <= 0) {
      tailp += lastp;
      tie_ct += (cmp_result == 0);
    }
  }
  // Sum away from center, until sums stop changing.
  while (1) {
    m11 += 1;
    m22 += 1;
    lastp_tail *= m12 * m21 / (m11 * m22);
    const double preaddp = tailp;
    tailp += lastp_tail;
    if (tailp == preaddp) {
      break;
    }
    m12 -= 1;
    m21 -= 1;
  }
  if (midp) {
    tailp -= S_CAST(double, tie_ct) * 0.5;
  }
  *resultp = starting_lnprob + log(tailp);
  return 0;
}

// obs_m11 + obs_m12 + obs_m21 + obs_m22 assumed to be <2^31.
// Just returns 0 if the log-p-value > 1 - 2^{-54} since additional precision
// in that direction is expected to be irrelevant.  (Reverse the test direction
// when you do actually want that precision.)
double Fisher22OneSidedLnP(uint32_t obs_m11, uint32_t obs_m12, uint32_t obs_m21, uint32_t obs_m22, uint32_t m11_is_greater_alt, uint32_t midp) {
  // Normalize.
  if (obs_m11 > obs_m22) {
    swap_u32(&obs_m11, &obs_m22);
  }
  if (obs_m12 > obs_m21) {
    swap_u32(&obs_m12, &obs_m21);
  }
  // Flipping m11<->m12 and m21<->m22 also flips the direction of the
  // alternative hypothesis.  So we flip on m11-is-greater alternative
  // hypothesis here to allow the rest of the code to assume m11-is-less.
  if (m11_is_greater_alt) {
    swap_u32(&obs_m11, &obs_m12);
    swap_u32(&obs_m21, &obs_m22);
  }
  const double obs_m11d = u31tod(obs_m11);
  const double obs_m12d = u31tod(obs_m12);
  const double obs_m21d = u31tod(obs_m21);
  const double obs_m22d = u31tod(obs_m22);
  double m11 = obs_m11d;
  double m12 = obs_m12d;
  double m21 = obs_m21d;
  double m22 = obs_m22d;
  if (S_CAST(uint64_t, obs_m11) * obs_m22 >= S_CAST(uint64_t, obs_m12) * obs_m21) {
    // We're to the right of the mode.
    // Start by computing an upper bound on the right-sum, and then iterating
    // leftward until we either know the p-value > 1 - 2^{-54} (at which point
    // we return 0), or remaining left likelihoods are smaller than the
    // precision limit.
    const double first_right_ratio = m12 * m21 / ((m11 + 1) * (m22 + 1));
    // 1 + r + r^2 + ... = 1 / (1-r)
    // const double right_upper_bound = 1.0 / (1 - first_right_ratio);

    // Rescale our starting lastp so that we overflow to INFINITY when we'd
    // want to early-exit and return 0; this saves us a comparison in the loop.
    const double left_rescale = (DBL_MAX / (1LL << 54)) * (1 - first_right_ratio) ;
    double lastp = left_rescale;
    double left_sum = 0;
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
    m11 = obs_m11d + 1;
    m12 = obs_m12d - 1;
    m21 = obs_m21d - 1;
    m22 = obs_m22d + 1;
    lastp = first_right_ratio;
    double right_sum = 1 + first_right_ratio;
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
    return log(1 - (right_sum - 0.5 * u31tod(midp)) / (right_sum + left_sum));
  }
  // We're to the left of the mode, and are responsible for tiny p-values.
  // If we're close enough to the mode that a simple left_sum / (left_sum +
  // right_sum) calculation doesn't risk overflow with the initial
  // relative-likelihood set to 1, just do that.
  // Otherwise... if the problem instance isn't *that* large, we could use
  // Lfact() to compute the starting log-likelihood and eat a big catastrophic
  // cancellation error, but I'll start with just the slow-and-accurate
  // calculation.
  const double m22_minus_m11 = obs_m22d - obs_m11d;
  const double m1x = obs_m11d + obs_m12d;
  const double m2x = obs_m21d + obs_m22d;
  const double mx1 = obs_m11d + obs_m21d;
  const double mxx = m1x + m2x;
  const double modal_m11 = m1x * mx1 / mxx;
  if (modal_m11 - obs_m11 <= 172) {
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
    m11 = obs_m11d;
    m12 = obs_m12d;
    m21 = obs_m21d;
    m22 = obs_m22d;
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
    return log((left_sum - 0.5 * u31tod(midp)) / (left_sum + right_sum));
  }
  // TODO
  return 0;
}

#ifdef __cplusplus
}
#endif
