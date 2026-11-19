#ifndef __FISHER_H__
#define __FISHER_H__

// Fisher's Exact Test library, copyright (C) 2013-2026 Christopher Chang.
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

// The main ideas behind these Fisher's Exact Test and hypergeometric-function
// implementations are:
// 1. It is expensive to evaluate the (log-)probability of a contingency table
//    from scratch, but cheap to compute the likelihood-ratio between a pair of
//    adjacent contingency tables.  Both Fisher's exact test and the
//    {p,q}hyper() functions involve summing probabilities of many adjacent
//    contingency tables; we express these sums as
//      <probability of starting table> * <multiplier>
//    and calculate the multiplier by repeatedly applying the likelihood-ratio
//    formula.
// 2. Log-probabilities drop off at a faster-than-quadratic rate as one moves
//    away from the mode.  When relative errors around 2^{-53} (float64) or
//    even 2^{-24} (float32) are acceptable, it's ok to ignore contingency
//    tables with probabilities less than roughly that multiple of the starting
//    contingency table's probability; for common large cases, we can ignore
//    most tables.
// 3. For 2x3 and larger Fisher's Exact Tests, we also take advantage of
//    combinatorial identities involving sums of entire lines/planes/etc. of
//    contingency table probabilities.  For 2xk and 3xk tables, we search
//    outward from the modal table, so that only the subproblems containing
//    tables more probable than the observed table need to be visited.
// 4. With the help of the QD high-precision library, we can evaluate the
//    log-probability of a single contingency table to better-than-float64
//    precision, and accumulate partial sums with better-than-float64 precision
//    when worthwhile.  This lets us achieve much higher accuracy *and* much
//    higher speed than widely-used implementations at the same time!
// (Implementations of the Mehta/Patel network algorithm get similar mileage
// from ideas #1 and #3, but take less advantage of #2 and #4.)
// 5. For larger 2xk and 3xk tables, we also borrow the network algorithm's
//    main idea, merging partial tables which leave the same subproblem
//    behind, in two ways:
//    - Meet-in-the-middle: the columns are split into two groups; each
//      nonempty configuration of the outer group is recorded as a query
//      against the subproblem it leaves behind (determined by the inner row
//      totals), and all queries against the same subproblem are answered
//      with a single enumeration of its inner tables.  This reduces the
//      exponent of the search cost by roughly a factor of two, and works
//      well with large cell counts.
//    - Memoization of subproblem results, keyed by the exact factorial
//      product of the fixed cells.  This works well with small cell counts,
//      where such products coincide frequently.  Unlike FEXACT's
//      tolerance-based merging, it doesn't introduce any error.
//    Since it's hard to predict in advance which of these is better, the
//    memoized search is tried first, with a work budget proportional to the
//    predicted meet-in-the-middle cost.

#ifdef __cplusplus
namespace plink2 {
#endif

double Fisher22TwoSidedP(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, int32_t midp, uint32_t logp);

double Fisher22TwoSidedPEx(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, td_real odds_tdr, int32_t midp, uint32_t logp);

double Fisher22OddsRatio(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22);

double Fisher22OddsRatioQuantileMatch(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double target_p);

// (low_p, high_p) should be the target quantiles.
// On exit, (*low_resultp, *high_resultp) is set to the odds ratios
// corresponding to those quantiles.
HEADER_INLINE void Fisher22OddsRatioCI(int64_t obs_m11, int64_t obs_m12, int64_t obs_m21, int64_t obs_m22, double low, double high, double* low_resultp, double* high_resultp) {
  // R's CI.
  *low_resultp = 1.0 / Fisher22OddsRatioQuantileMatch(obs_m12, obs_m11, obs_m22, obs_m21, low);
  *high_resultp = Fisher22OddsRatioQuantileMatch(obs_m11, obs_m12, obs_m21, obs_m22, 1.0 - high);
}

double Fisher23LnP(int32_t obs_m11, int32_t obs_m12, int32_t obs_m13, int32_t obs_m21, int32_t obs_m22, int32_t obs_m23, uint32_t midp);

// obs[] is in row-major order: obs[0..(col_ct-1)] is the first row, and
// obs[col_ct..(2*col_ct-1)] is the second row.  Entries must be nonnegative
// and sum to less than 2^31.
// Returns log-[mid]p-value, or NaN on allocation failure.
// If there are more than 3 nonempty columns, and a cheap upper bound on the
// log-p-value is less than ln_p_cutoff, that upper bound is returned instead.
// (The bound is ln(observed table probability) + ln(bound on number of tables
// with the same margins), and since the p-value is at least the observed
// table probability, it overshoots by less than the latter term, which is
// roughly (k-1) * ln(n).  So this mainly helps when the p-value is tiny.)
// Pass -DBL_MAX to always compute the exact value.
// With large cell counts, search cost is roughly O(w^{k/2}), where k is the
// number of nonempty columns and w is proportional to the standard deviation
// of a column's row-1 entry; with small cell counts, it's usually much
// lower.
// Memory use grows with the difficulty of the problem: it's usually well
// under 50 MiB, and internal caps bound it to roughly 250 MiB for tables
// with up to ~20 columns.  The caps are: 72 MiB for the memoization table;
// ~64 MiB plus ~1 MiB per column (2 MiB for 3xk) for a 2^19-query
// meet-in-the-middle batch; and 64 MiB and 32 MiB for the log-factorial and
// factorial-hash tables, which are proportional to the table total until
// then.  Problems approaching the bound take many seconds; e.g. a 3x5 table
// with total 450 and p ~= 0.02 takes ~15 s and ~180 MiB.
double Fisher2xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff);

// obs[] is a 3xk table in row-major order: obs[0..(col_ct-1)] is the first
// row, etc.  Entries must be nonnegative and sum to less than 2^31.  Other
// semantics are the same as for Fisher2xkLnP().  (Search cost is roughly
// O(w^k) for large cell counts, with w defined as for Fisher2xkLnP().  Memory
// use is as described for Fisher2xkLnP().)
double Fisher3xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff);

// obs[] is a 3x3 table in row-major order.  Equivalent to
// Fisher3xkLnP(obs, 3, midp, ln_p_cutoff).
double Fisher33LnP(const int32_t* obs, uint32_t midp, double ln_p_cutoff);

// Probable todos:
// - Larger jxk tables.  The meet-in-the-middle and memoization machinery is
//   mostly table-shape-agnostic, but the outward-from-the-mode search
//   currently relies on 2- and 3-row-specific mode updates.

// Practical limits of Fisher2xkLnP() and Fisher3xkLnP() (with ln_p_cutoff =
// -DBL_MAX), as of Nov 2026 on one core of an Apple M2 Pro.  Each entry is the
// approximate table total which can typically be handled within one second /
// within one minute.  Tables had moderately unequal margins, plus a deviation
// from independence scaled to produce the given p-value; times are medians
// over 3 tables.  Individual tables vary by a factor of several (much more
// with many columns, where the totals below are rough), and runtime grows
// quickly as the p-value shrinks further.
//
//   shape     p ~= 0.5          p ~= 1e-3         p ~= 1e-10
//   2x4       1.0M / 16M        320k / 5M         170k / 2.7M
//   2x6       37k  / 230k       9k   / 53k        3.7k / 21k
//   2x8       4.1k / 14k        1.1k / 3.8k       460  / 1.5k
//   2x10      1.1k / 2.7k       320  / 830        170  / 360
//   2x15      220  / 380        125  / 180        -
//   3x3       38k  / 310k       11k  / 89k        5k   / 43k
//   3x4       3k   / 14k        770  / 3.5k       290  / 1.4k
//   3x5       400  / 1.6k       120  / 420        110  / 190
//   3x6       210  / 730        95   / 240        -    / 140
//   3x8       110  / 170        75   / 95         -
//
// ("-" = such a small p-value isn't attainable at totals small enough.)
// 2x3 tables, handled by Fisher23LnP(), are fast all the way to the 2^31
// input limit: ~0.5 s at a total of 3e8 with p ~= 0.5, and ~13 s at 2e9 with
// p ~= 1e-10.  When the caller only needs to know whether the p-value is
// below some threshold, or only needs a non-log p-value (which underflows to
// zero below ~1e-323), passing an appropriate ln_p_cutoff can avoid the
// worst cases.

#ifdef __cplusplus
}
#endif

#endif  // __FISHER_H__
