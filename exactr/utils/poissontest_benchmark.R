#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))

description <- "This collects one-sample two-sided benchmark results for exactr::poisson.test()
and stats::poisson.test().  The latter takes O(r*T) time and memory, so it is
skipped for r*T >= 2^max_stats_exp."

poissontest_benchmark <- function(z, pow2s, num_trials_per_pow2, max_stats_exp) {
  warmup <- exactr::poisson.test(1, r = 2.5)
  warmup <- stats::poisson.test(1, r = 2.5)
  for (pow2 in pow2s) {
    m <- 2^pow2 - 0.5
    k <- round(m + z * sqrt(m))
    if (k < 0) {
      k <- 0
    }
    secs_base <- seconds_per_iter(function() exactr::poisson.test(k, r = m), num_trials_per_pow2)
    print_str <- paste0("r*T=(2^", pow2, ")-0.5: base=", format_g3(secs_base))
    if (pow2 < max_stats_exp) {
      secs_stats <- seconds_per_iter(function() stats::poisson.test(k, r = m), num_trials_per_pow2)
      print_str <- paste0(print_str, "  stats=", format_g3(secs_stats))
    }
    cat(print_str, " sec/iter\n", sep = "")
  }
}

main <- function() {
  cmd_args <- parse_commandline_args(description, list(
    list(short = "z", long = "z-score", type = "numeric", default = 0,
         help = "Event-count z-score to test."),
    list(short = "e", long = "exps", type = "string", default = "5,20,35",
         help = "Test r*T=2**<these values> - 0.5."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "m", long = "max-stats-exp", type = "integer", default = 25L,
         help = "Skip stats::poisson.test() for r*T >= 2**<this value>.")))
  poissontest_benchmark(cmd_args$z_score,
                        parse_range_string(cmd_args$exps),
                        cmd_args$number,
                        cmd_args$max_stats_exp)
}

main()
