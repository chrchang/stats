#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))
source(file.path(.script_dir, "mpfr_impls.R"))

description <- "Primary mode checks exactr::pbinom(approx=FALSE) and stats::pbinom()'s
accuracy against a simple MPFR-based implementation (slow!).

\"--bits 0\" mode checks exactr::pbinom(approx=TRUE) and stats::pbinom(),
treating exactr::pbinom(approx=FALSE) as ground truth."

pbinom_accuracy_test <- function(p, z, pow2s, num_trials_per_pow2, bits, logp) {
  pq <- p * (1 - p)
  for (pow2 in pow2s) {
    relerr_noapprox_ssq <- 0
    relerr_approx_ssq <- 0
    relerr_stats_ssq <- 0
    ns <- sample_pow2_tier(pow2, num_trials_per_pow2)
    for (n in ns) {
      stdev <- sqrt(n * pq)
      k <- round(n * p + z * stdev + 0.5)
      if (k < 0) {
        k <- 0
      } else if (k > n) {
        k <- n
      }
      got_noapprox <- flush_if_denormal(exactr::pbinom(k, n, p, log.p = logp))
      got_approx <- flush_if_denormal(exactr::pbinom(k, n, p, log.p = logp, approx = TRUE))
      got_stats <- stats::pbinom(k, n, p, log.p = logp)
      if (bits == 0) {
        want <- got_noapprox
      } else {
        want <- flush_if_denormal(binom_cdf(k, n, p, bits, logp))
        relerr <- compute_relerr(got_noapprox, want)
        relerr_noapprox_ssq <- relerr_noapprox_ssq + relerr * relerr
      }
      relerr <- compute_relerr(got_approx, want)
      relerr_approx_ssq <- relerr_approx_ssq + relerr * relerr
      relerr <- compute_relerr(got_stats, want)
      relerr_stats_ssq <- relerr_stats_ssq + relerr * relerr
    }
    num_trials <- length(ns)
    approx_rms <- sqrt(relerr_approx_ssq / num_trials)
    stats_rms <- sqrt(relerr_stats_ssq / num_trials)
    print_str <- paste0("n in [2^", pow2, ", 2^", pow2 + 1, "): ")
    if (bits != 0) {
      noapprox_rms <- sqrt(relerr_noapprox_ssq / num_trials)
      print_str <- paste0(print_str, "errRMS=", format_g3(noapprox_rms), "  ")
    }
    cat(print_str, "approxErrRMS=", format_g3(approx_rms), "  statsErrRMS=", format_g3(stats_rms), "\n", sep = "")
  }
}

main <- function() {
  cmd_args <- parse_commandline_args(description, list(
    list(short = "p", long = "succ-prob", type = "numeric", default = 0.5,
         help = "Binomial distribution success-probability to test."),
    list(short = "z", long = "z-score", type = "numeric", default = 0,
         help = "Success-count z-score to test."),
    list(short = "e", long = "exps", type = "string", default = "5,20,35",
         help = "Test n~=2**<these values>."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "b", long = "bits", type = "integer", default = 256L,
         help = "MPFR precision; or if 0, treat approx=FALSE as ground truth."),
    list(short = "l", long = "logp", type = "flag",
         help = "Test log.p=TRUE."),
    list(short = "s", long = "seed", type = "integer", default = 1L,
         help = "RNG seed.")))
  set.seed(cmd_args$seed)
  pbinom_accuracy_test(cmd_args$succ_prob,
                       cmd_args$z_score,
                       parse_range_string(cmd_args$exps),
                       cmd_args$number,
                       cmd_args$bits,
                       cmd_args$logp)
}

main()
