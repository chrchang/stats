# exactr 0.1.3 (development)
* fisher.test() on 3xk and kx3 tables now usually uses a network algorithm
  first (merging partial tables that leave the same subproblem, with bounds on
  each subproblem's most and least probable completions), which is much faster
  on most tables taking more than a few milliseconds, especially wider ones,
  with the same exact tie handling, midp support, and accuracy.  Its memory
  limit is max(4 * workspace, 256 * 2^20) bytes (256 MiB by default); if the
  limit would be exceeded, the previous algorithm is used.
* fisher.test() now also handles tables with at least 4 rows and 4 columns
  (at most 8 nonempty rows or columns) itself, with the same network
  algorithm, instead of forwarding them to stats::fisher.test(); this is
  usually faster, supports midp = TRUE, resolves ties exactly instead of with
  a ~1e-7 relative tolerance, and handles some tables where FEXACT fails (e.g.
  with a hash key overflow).  If the memory limit would be exceeded, or the
  table is larger, it is still forwarded.
* fisher.test() on 2xk and kx2 tables with at least 8 nonempty columns and
  small counts (total at most 288, and at most 25 times the number of
  columns) also tries the network algorithm first, with the same memory
  limit; it is usually 5-100x faster on such tables, and a larger workspace can help a lot with moderate
  counts.

# exactr 0.1.2 (2026-11-18)
* New poisson.test().
* fisher.test() now handles all 2xk, kx2, 3xk, and kx3 tables itself (when
  simulate.p.value = FALSE), instead of forwarding 2xk tables with more than
  6 nonempty columns and 3xk tables with k > 3 to stats::fisher.test().  The
  new algorithm combines meet-in-the-middle and exact memoization ideas from
  the network algorithm with the previous approach; it's usually much faster
  than stats::fisher.test() on these tables (often by 10-100x with larger
  counts, and avoiding its workspace errors), manages its own memory (usually
  well under 50 MB, at most ~250 MB with up to ~20 columns) instead of
  requiring a 'workspace' setting, supports midp = TRUE, and computes p-values
  to nearly full double precision.
* fisher.test()'s default workspace (used only for tables handled by
  stats::fisher.test()'s FEXACT, i.e. those with at least 4 rows and 4
  columns) is now 64000000 (256 MB) instead of 200000, in line with the memory
  limit of the new 2xk/3xk algorithm.
* fisher.test() on 2xk tables with 4+ nonempty columns and 3x3 tables is more
  accurate with larger counts (previously, relative error could reach ~1e-11
  at n = 1000).
* binom.test(alternative="greater", midp=TRUE) now returns
  P[X > x] + 0.5 * P[X = x]; it previously returned
  P[X >= x-1] - 0.5 * P[X = x-1].
* binom.test(alternative="greater") now returns 1 when x = 0 (previously 0).
* binom.test(midp=TRUE) now reports method "Exact binomial test (midp)".

# exactr 0.1.0 (2026-10-28)
Initial CRAN submission.
