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
//    when worthwhile.
// 5. For larger 2xk and 3xk tables, the outward search also borrows the
//    network algorithm's main idea, merging partial tables which leave the
//    same subproblem behind, in two ways:
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
// 6. Most 3xk tables, 2xk tables with many columns and small counts, and all
//    tables with 4+ rows and columns are handled by a full network algorithm
//    (fisher_network.cc), which combines it with ideas #1, #2, and #4: nodes
//    are sorted vectors of remaining row totals; each node's list of "past
//    values" is merged on exact equality of factorial products instead of
//    within a tolerance; exact transportation modes, plus cheaper bounds,
//    let whole blocks of a sorted list be counted in closed form or dropped;
//    the last two columns are enumerated outward from the mode; and all
//    weights come from double-double log-factorials.  It is usually much
//    faster than the outward search on the tables it is used for, and than
//    FEXACT on 4x4-6x6 tables (see the notes after FisherRxcLnPEx(), and the
//    practical-limits tables at the end of this file).
//
// Current status (Nov 2026):
// - Two-sided tests, where a table counts as at least as extreme when its
//   probability is <= the observed table's (one-sided 2x2 tests are
//   hypergeometric tail probabilities; see hypergeom.h); log-p-values and
//   mid-p variants are available.  Ties are resolved exactly (FEXACT, as used
//   by R's fisher.test(), counts tables within a relative tolerance of ~3e-7
//   of the observed probability as ties), and p-values are accurate to
//   nearly float64 precision (stock FEXACT's running sum can be off by
//   ~1e-8 relative on large tables).
// - 2x2: Fisher22TwoSidedP[Ex](), with the conditional MLE odds ratio and
//   its confidence interval.  2x3: Fisher23LnP(), fast up to the 2^31 input
//   limit.
// - 2xk and 3xk: Fisher2xkLnP[Ex]() and Fisher3xkLnP[Ex](), which choose
//   between the outward search and the network algorithm and fall back from
//   the latter automatically, so they always return an exact p-value given
//   enough time.
// - Larger r x c tables, with at most kFisherRxcMaxMinDim (8) nonempty rows
//   or columns: FisherRxcLnPEx().  With 4+ nonempty rows and columns it uses
//   the network algorithm alone, so it can fail with kFisherRxcMemLimit.
// - The *Ex() functions take a network memory limit and a cancellation hook.
//   The network algorithm keeps all of its state per call, so concurrent
//   calls are safe.

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
// Tables with at least 8 nonempty columns and small counts (total at most
// 288, and at most 25 times the number of nonempty columns), unless the
// observed table is close to the mode, are first handled by the network
// algorithm described for Fisher3xkLnP() (with memory limit
// kFisherNetDefaultMemLimit here); it is usually 5-100x faster on such tables,
// and a failed attempt costs at most ~1 s at the default limit.  If its memory
// limit would be exceeded, the search described above is used.  (With 2-row
// tables the network's memory need grows quickly with the counts per column:
// in Nov 2026 tests, many 2x10-2x32 tables with totals ~200-300 needed limits
// of 256-512 MiB and took under 1 s, and a few needed ~1 GiB and took ~3.5 s,
// vs. 7 s to many minutes for the search; tables with totals ~400+ mostly
// needed more than 1 GiB.)
double Fisher2xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff);

// Default memory limit for the network algorithm used by Fisher2xkLnPEx(),
// Fisher3xkLnPEx(), and FisherRxcLnPEx() (256 MiB).
static const uint64_t kFisherNetDefaultMemLimit = 268435456;

// Cancellation hook for the *Ex() functions below.  Runtime is hard to
// predict in advance (cheap cost estimates are often off by 10x or more, and
// the memoized search frequently finishes in a tiny fraction of the
// estimate), so callers wanting a time limit or interruptibility should poll
// instead.  If check_fn is non-null, it's called with check_arg periodically
// during the search (roughly every 1-10 ms on a modern CPU); if it returns
// nonzero, the search is abandoned and NaN is returned.  (NaN is also
// returned on allocation failure; the caller can distinguish the two cases by
// tracking whether its check function returned nonzero.)  check_fn must not
// longjmp out (e.g. via R_CheckUserInterrupt()), since that leaks memory.
typedef uint32_t (*FisherCheckFn)(void* check_arg);

// Fisher2xkLnP() with a network memory limit (see Fisher3xkLnPEx(); 0 means
// no limit, which is dangerous: see the warning there) and a cancellation
// hook.
double Fisher2xkLnPEx(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff, uint64_t mem_limit_bytes, FisherCheckFn check_fn, void* check_arg);


// obs[] is a 3xk table in row-major order: obs[0..(col_ct-1)] is the first
// row, etc.  Entries must be nonnegative and sum to less than 2^31.  Other
// semantics are the same as for Fisher2xkLnP().
// Unless the observed table is close to the mode (where the outward search
// below is cheap), tables with 3 or more nonempty columns are first handled by
// a network algorithm (fisher_network.cc), which merges partial tables leaving
// the same subproblem and bounds each subproblem's most and least probable
// completions; this is usually much faster than the outward search, especially
// for wider tables (e.g. 1.2-1.5x on 3x3 and 1.1-11x on 3x4 tables taking more
// than a few milliseconds, and up to several hundred times on 3x6-3x8 tables
// with small p-values).  If the network's workspace would exceed the memory
// limit (kFisherNetDefaultMemLimit for this function), or the table total is
// at least 2^20, the outward search with meet-in-the-middle and
// memoization (cost roughly O(w^k) for large cell counts, with w defined as
// for Fisher2xkLnP(); memory as described for Fisher2xkLnP()) is used instead.
double Fisher3xkLnP(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff);

// Fisher3xkLnP() with a network memory limit and a cancellation hook.
// mem_limit_bytes bounds the network algorithm's main workspace: its
// past-value lists, run descriptors, node records, hash maps, merge/sort
// buffers, and log-factorial tables; 0 means no limit.  Its large buffers are
// mapped from (and returned to) the OS directly, so the allocator does not
// keep freed blocks around; small allocations and the process's other memory
// are not counted.  In Nov 2026 tests, peak resident memory was within ~5%
// (at most ~30 MiB) of the limit with 200 MiB to 1 GiB limits.
// Under memory pressure the network first merges pending partial tables early,
// then switches some subproblems to a slower list-free mode, and finally falls
// back to the outward search, whose memory is bounded separately (see
// Fisher2xkLnP()).
// Warning: with no limit (or one above the free physical memory), a hard
// table makes the network allocate until the system runs out of memory; in
// tests, wide tables with moderate counts consumed several GiB within ~10
// seconds, well before typical time limits, and often before a user would
// notice and cancel.
// Depending on the OS, that can make the whole machine unresponsive, crash
// it, or get this or other processes killed, instead of failing cleanly.  Use
// a finite limit comfortably below the available physical memory.
double Fisher3xkLnPEx(const int32_t* obs, uint32_t col_ct, uint32_t midp, double ln_p_cutoff, uint64_t mem_limit_bytes, FisherCheckFn check_fn, void* check_arg);

// obs[] is a 3x3 table in row-major order.  Equivalent to
// Fisher3xkLnP(obs, 3, midp, ln_p_cutoff).
double Fisher33LnP(const int32_t* obs, uint32_t midp, double ln_p_cutoff);

// Status codes for FisherRxcLnPEx().
enum {
  kFisherRxcOk = 0,
  // The network algorithm's workspace would exceed mem_limit_bytes.
  kFisherRxcMemLimit = 1,
  // More than kFisherRxcMaxMinDim nonempty rows and columns, or (very rarely)
  // a case the network algorithm can't resolve.
  kFisherRxcUnsupported = 2,
  // check_fn returned nonzero, or an allocation failed.
  kFisherRxcFailed = 3
};

// Largest supported min(nonempty rows, nonempty columns) for
// FisherRxcLnPEx().  This is a compile-time cap, not a correctness limit: the
// network algorithm's per-node and scratch arrays are sized for it
// (kFisherNetMaxRows in fisher_network.cc).  Raising it is mostly mechanical,
// but tables with many rows *and* columns quickly become intractable for any
// exact network method: the number of subproblems grows roughly like
// R^(r-1) and the number of column configurations per subproblem like
// c^(r-1) / (r-1)!, and the exact group minimum used for pruning checks
// j * 2^(j-1) vertices for j free rows (1024 at 8 rows, ~500k at 16), so
// only small, sparse 9+ x 9+ tables would benefit.  (FEXACT usually fails on
// such tables as well.)
static const uint32_t kFisherRxcMaxMinDim = 8;

// obs[] is a row_ct x col_ct table in row-major order.  Entries must be
// nonnegative and sum to less than 2^31.  Returns the log-[mid]p-value of the
// two-sided Fisher's exact test (tables count as at least as extreme when
// their probability is <= the observed table's; ties are resolved exactly).
// Empty rows and columns are removed first.  Tables with at most 3 nonempty
// rows or columns are handled by Fisher22TwoSidedP(), Fisher2xkLnPEx(), or
// Fisher3xkLnPEx(); larger ones (with at most kFisherRxcMaxMinDim nonempty
// rows or columns) by the network algorithm in fisher_network.cc, with
// mem_limit_bytes as described for Fisher3xkLnPEx() but no fallback.
// ln_p_cutoff and check_fn work as for Fisher2xkLnPEx().  On failure, NaN is
// returned and *status_out says why.  (For tables whose p-value is very close
// to 1, the network sums the p-value directly, so the log-p-value has absolute
// error around 1e-15 rather than full relative precision.)
double FisherRxcLnPEx(const int32_t* obs, uint32_t row_ct, uint32_t col_ct, uint32_t midp, double ln_p_cutoff, uint64_t mem_limit_bytes, FisherCheckFn check_fn, void* check_arg, uint32_t* status_out);

// Possible todos:
// - The network's double-double low parts cost ~25-30% memory; storing them
//   only where needed would let it fit somewhat larger tables.

// Practical limits of Fisher2xkLnP(), Fisher3xkLnP(), and (for 4+ rows)
// FisherRxcLnPEx() with mem_limit_bytes = kFisherNetDefaultMemLimit (all with
// ln_p_cutoff = -DBL_MAX), as of Nov 2026 on an Apple M2 Pro (one core per
// table).  Each entry is the approximate table total which can typically be
// handled within one second / within one minute.  Tables had moderately
// unequal margins, plus a deviation from independence scaled to produce the
// given p-value; times are medians over 3 tables.  Individual tables vary by
// a factor of several (much more with many columns, where the totals below
// are rough), and runtime grows quickly as the p-value shrinks further.
//
//   shape     p ~= 0.5          p ~= 1e-3         p ~= 1e-10
//   2x4       940k / 15M        290k / 4.7M       150k / 2.6M
//   2x6       57k  / 280k       13k  / 61k        4.7k / 23k
//   2x8       4.6k / 15k        1.2k / 4.0k       480  / 1.6k
//   2x10      1.1k / 2.9k       350  / 860        200  / 390
//   2x15      280  / 400        240  / 270        -
//   2x20      280  / 320        240  / 270        -
//   3x3       42k  / 350k       11k  / 98k        5.1k / 46k
//   3x4       5.5k / 16k        1.4k / 3.9k       690  / 1.5k
//   3x5       1.4k / 3.5k       440  / 1.2k       220  / 550
//   3x6       520  / 750        220  / 270        160  / 180
//   3x8       260  / 290        160  / 180        -
//   3x10      220  / 270        150  / 170        -
//   4x4       470  / 470m       200m              -
//   4x5       270m              130  / 170m       -
//   4x6       170  / 240m       96   / 150m       -
//   4x8       150  / 180m       -                 -
//   5x5       130  / 160m       -                 -
//   5x6       115  / 145m       -                 -
//   6x6       86   / 105m       -                 -
//
// ("-" = such a small p-value isn't attainable at totals small enough.  "m" =
// limited by memory rather than time: above roughly this total, the network
// algorithm usually exceeds the 256 MiB default limit before the time is up,
// and FisherRxcLnPEx() returns kFisherRxcMemLimit; a single "m" entry means
// this happens in under a second.  A larger mem_limit_bytes extends these
// shapes somewhat, at the cost of the memory.  2xk and 3xk tables fall back
// to a search with separately bounded memory instead.)
// 2x3 tables, handled by Fisher23LnP(), are fast all the way to the 2^31
// input limit: ~0.5 s at a total of 3e8 with p ~= 0.5, and ~13 s at 2e9 with
// p ~= 1e-10.  When the caller only needs to know whether the p-value is
// below some threshold, or only needs a non-log p-value (which underflows to
// zero below ~1e-323), passing an appropriate ln_p_cutoff can avoid the
// worst cases.

// For comparison, the same measurements for stock FEXACT, as called by R
// 4.7.0's fisher.test() (r90935) with workspace = 2^26 (256 MiB, vs. the
// 800 KB default) and the default mult = 30:
//
//   shape     p ~= 0.5          p ~= 1e-3         p ~= 1e-10
//   2x4       56k  / 88kx       22k  / 92kx       13k  / 81kx
//   2x6       2.4k / 8.8km      780  / 2.0km      400  / 760m
//   2x8       630  / 970m       250  / 310m       160  / 200m
//   2x10      310  / 390m       200  / 240m       -    / 190m
//   2x15      250  / 360m       170  / 240m       -
//   2x20      100x              -                 -
//   3x3       1.4k / 3.7kx      960  / 3.4k       750  / 2.6k
//   3x4       520  / 1.3k       260  / 600        200  / 370
//   3x5       270  / 480        140  / 260        -    / 160m
//   3x6       190  / 300m       110  / 180m       -
//   3x8       160  / 250m       -    / 160m       -
//   3x10      140  / 180x       -                 -
//   4x4       140  / 240        96   / 160        -
//   4x5       120  / 190        85   / 130        -
//   4x6       110  / 160        -                 -
//   4x8       94x               -                 -
//   5x5       84   / 140        -                 -
//   5x6       80   / 120        -                 -
//   6x6       66   / 70x        -                 -
//
// ("m" = above roughly this total, most tables exhaust the workspace (LDSTP
// or LDKEY too small) before the time is up; "x" = most fail with an integer
// overflow ("error 501", a hash key above INT_MAX, or "integer overflow in
// exact computation"), which a larger workspace doesn't fix; a single
// "m"/"x" entry means this happens in under a second.  "-" as above, also
// counting totals where FEXACT already fails.)  So with the same memory,
// the functions here handle larger totals for every shape: ~1.1-2x on 2x15,
// 3x8, 3x10, and 4x5 through 6x6 (~3x on 2x20), 2-15x on 2x8, 2x10, 3x4-3x6,
// and 4x4, and 15-170x on 2x4, 2x6, and 3x3.  FEXACT is also less accurate: it
// counts tables within a relative tolerance of ~3e-7 of the observed
// probability as ties, and its running sum can be off by ~1e-8 relative on
// large tables.

#ifdef __cplusplus
}
#endif

#endif  // __FISHER_H__
