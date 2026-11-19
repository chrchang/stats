#ifndef __POIS_H__
#define __POIS_H__

// Poisson Exact Test library, copyright (C) 2026 Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_base.h"
#include "plink2_highprec.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// This follows the same strategy as binom.h: evaluate the log-probability of
// a single outcome to better-than-float64 precision (with Loader's
// algorithm), and multiply it by a tail-sum of relative likelihoods.  The
// tail-sum is computed by
// - directly summing adjacent-outcome likelihood ratios (lambda/(k+1) or
//   k/lambda) when that converges quickly,
// - the n -> infinity limit of the Aroian/DiDonato/Morris continued fraction
//   used by binom_tail_lik_bfrac() (upper tail), or Legendre's continued
//   fraction for the incomplete gamma function (lower tail; this is the other
//   n -> infinity limit of the same continued fraction), when we're more than
//   a few standard deviations from the mean, and
// - the n -> infinity limit of BASYM (i.e. Temme's uniform asymptotic
//   expansion of the incomplete gamma function ratio) near the mean.
// Overall accuracy is similar to that of PbinomApprox() and BinomTwoSidedP().
//
// Unlike the binomial case, the Poisson distribution is not symmetric, so the
// two-sided test is not normalized to a single tail.
//
// lambda is passed as a dd_real, since R poisson.test() computes it as the
// product (r * T) of two doubles, which a dd_real represents exactly.

// Assumes 0 <= k < 2^53 integer, 0 < lambda < 2^960.
double PoisMass(double k, dd_real lambda_ddr, uint32_t logp);

// Returns P[X <= obs_k] if complement is false, P[X > obs_k] otherwise.  If
// midp is true, 0.5 * P[X = obs_k] is subtracted from the former / added to
// the latter.
// Assumes 0 <= lambda < 2^960 (lambda=0 ok) and obs_k is an integer.
// obs_k >= 2^53 is only supported when lambda < 2^51.
double PpoisApprox(double obs_k, dd_real lambda_ddr, int32_t complement, int32_t midp, uint32_t logp);

// High-accuracy variant of PpoisApprox() (without midp support).  Should
// consistently achieve <1 ULP relative error, like Pbinom().  Same domain as
// PpoisApprox().
double Ppois(double obs_k, dd_real lambda_ddr, int32_t complement, uint32_t logp);

// One-sided exact test p-value; greater_alt indicates alternative="greater".
HEADER_INLINE double PoisOneSidedP(double obs_k, dd_real lambda_ddr, uint32_t greater_alt, int32_t midp, uint32_t logp) {
  const int32_t k_decr = greater_alt && (!midp);
  if (k_decr && (obs_k == 0)) {
    return logp? 0.0 : 1.0;
  }
  return PpoisApprox(obs_k - k_decr, lambda_ddr, greater_alt, midp, logp);
}

// Two-sided ("minlike") exact test p-value: sum of probabilities of all
// outcomes with probability <= P[X = obs_k].  (Ties are resolved exactly,
// rather than with R poisson.test()'s 1 + 1e-7 relative tolerance.)  When midp
// is true, outcomes with probability exactly equal to P[X = obs_k] are given
// half weight.
// Assumes 0 <= obs_k < 2^52, 0 <= lambda < 2^50.
double PoisTwoSidedP(int64_t obs_k, dd_real lambda_ddr, int32_t midp, uint32_t logp);

#ifdef __cplusplus
}
#endif

#endif  // __POIS_H__
