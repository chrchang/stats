#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))

description <- "This collects two-sided benchmark results for exactr::binom.test() and
stats::binom.test().  The latter allocates O(n) memory, so it is skipped for
n >= 2^max_stats_exp."

binomtest_benchmark <- function(p, z, pow2s, num_trials_per_pow2, max_stats_exp) {
  pq <- p * (1 - p)
  warmup <- exactr::binom.test(1, 2, 0.5)
  warmup <- stats::binom.test(1, 2, 0.5)
  for (pow2 in pow2s) {
    n <- 2^pow2 - 1
    stdev <- sqrt(n * pq)
    k <- round(n * p + z * stdev)
    if (k < 0) {
      k <- 0
    } else if (k > n) {
      k <- n
    }
    secs_base <- seconds_per_iter(function() exactr::binom.test(k, n, p), num_trials_per_pow2)
    print_str <- paste0("n=(2^", pow2, ")-1: base=", format_g3(secs_base))
    if (pow2 < max_stats_exp) {
      secs_stats <- seconds_per_iter(function() stats::binom.test(k, n, p), num_trials_per_pow2)
      print_str <- paste0(print_str, "  stats=", format_g3(secs_stats))
    }
    cat(print_str, " sec/iter\n", sep = "")
  }
}

main <- function() {
  cmd_args <- parse_commandline_args(description, list(
    list(short = "p", long = "succ-prob", type = "numeric", default = 0.5,
         help = "Binomial distribution success-probability to test."),
    list(short = "z", long = "z-score", type = "numeric", default = 0,
         help = "Success-count z-score to test."),
    list(short = "e", long = "exps", type = "string", default = "5,20,35",
         help = "Test n=2**<these values> - 1."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "m", long = "max-stats-exp", type = "integer", default = 27L,
         help = "Skip stats::binom.test() for n >= 2**<this value>.")))
  binomtest_benchmark(cmd_args$succ_prob,
                      cmd_args$z_score,
                      parse_range_string(cmd_args$exps),
                      cmd_args$number,
                      cmd_args$max_stats_exp)
}

main()
