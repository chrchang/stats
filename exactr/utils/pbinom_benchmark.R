#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))

description <- "This collects benchmark results for exactr::pbinom() and stats::pbinom()."

pbinom_benchmark <- function(p, z, pow2s, num_trials_per_pow2, logp) {
  pq <- p * (1 - p)
  warmup <- exactr::pbinom(1, 2, 0.5)
  warmup <- stats::pbinom(1, 2, 0.5)
  for (pow2 in pow2s) {
    n <- 2^pow2 - 1
    stdev <- sqrt(n * pq)
    k <- round(n * p + z * stdev)
    if (k < 0) {
      k <- 0
    } else if (k > n) {
      k <- n
    }
    secs_noapprox <- seconds_per_iter(function() exactr::pbinom(k, n, p, log.p = logp), num_trials_per_pow2)
    secs_approx <- seconds_per_iter(function() exactr::pbinom(k, n, p, log.p = logp, approx = TRUE), num_trials_per_pow2)
    secs_stats <- seconds_per_iter(function() stats::pbinom(k, n, p, log.p = logp), num_trials_per_pow2)
    cat("n=(2^", pow2, ")-1: base=", format_g3(secs_noapprox), "  approx=", format_g3(secs_approx), "  stats=", format_g3(secs_stats), " sec/iter\n", sep = "")
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
    list(short = "l", long = "logp", type = "flag",
         help = "Test log.p=TRUE.")))
  pbinom_benchmark(cmd_args$succ_prob,
                   cmd_args$z_score,
                   parse_range_string(cmd_args$exps),
                   cmd_args$number,
                   cmd_args$logp)
}

main()
