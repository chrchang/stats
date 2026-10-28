#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))
source(file.path(.script_dir, "mpfr_impls.R"))

description <- "Primary mode checks exactr:::phyper(approx=FALSE) and stats::phyper()'s
accuracy against a simple MPFR-based implementation (slow!).

\"--bits 0\" mode checks exactr:::phyper(approx=TRUE) and stats::phyper(),
treating exactr:::phyper(approx=FALSE) as ground truth."

phyper_accuracy_test <- function(ab, ac, z, pow2s, num_trials_per_pow2, bits, logp, omit_stats) {
  for (pow2 in pow2s) {
    ns <- sample_pow2_tier(pow2, num_trials_per_pow2)
    relerr_noapprox_ssq <- 0
    relerr_approx_ssq <- 0
    relerr_stats_ssq <- 0
    for (n in ns) {
      tbl <- make_22_table(n, ab, ac, z)
      a <- tbl[1]
      cur_ab <- tbl[1] + tbl[2]
      cur_ac <- tbl[1] + tbl[3]
      got_noapprox <- flush_if_denormal(exactr:::phyper(a, cur_ac, n - cur_ac, cur_ab, log.p = logp))
      got_approx <- flush_if_denormal(exactr:::phyper(a, cur_ac, n - cur_ac, cur_ab, log.p = logp, approx = TRUE))
      if (bits == 0) {
        want <- got_noapprox
      } else {
        want <- flush_if_denormal(hypergeom_cdf(a, n, cur_ab, cur_ac, bits, logp))
        relerr <- compute_relerr(got_noapprox, want)
        relerr_noapprox_ssq <- relerr_noapprox_ssq + relerr * relerr
      }
      relerr <- compute_relerr(got_approx, want)
      relerr_approx_ssq <- relerr_approx_ssq + relerr * relerr
      if (!omit_stats) {
        got_stats <- stats::phyper(a, cur_ac, n - cur_ac, cur_ab, log.p = logp)
        relerr <- compute_relerr(got_stats, want)
        relerr_stats_ssq <- relerr_stats_ssq + relerr * relerr
      }
    }
    num_trials <- length(ns)
    print_str <- paste0("n in [2^", pow2, ", 2^", pow2 + 1, "): ")
    if (bits > 0) {
      noapprox_rms <- sqrt(relerr_noapprox_ssq / num_trials)
      print_str <- paste0(print_str, "errRMS=", format_g3(noapprox_rms), "  ")
    }
    approx_rms <- sqrt(relerr_approx_ssq / num_trials)
    print_str <- paste0(print_str, "approxErrRMS=", format_g3(approx_rms))
    if (!omit_stats) {
      stats_rms <- sqrt(relerr_stats_ssq / num_trials)
      print_str <- paste0(print_str, "  statsErrRMS=", format_g3(stats_rms))
    }
    cat(print_str, "\n", sep = "")
  }
}

main <- function() {
  cmd_args <- parse_commandline_args(description, list(
    list(short = "r", long = "row-prop", type = "numeric", default = 0.5,
         help = "Proportion of total in first row."),
    list(short = "c", long = "col-prop", type = "numeric", default = 0.5,
         help = "Proportion of total in first column."),
    list(short = "z", long = "z-score", type = "numeric", default = 0,
         help = "Z-score to test."),
    list(short = "e", long = "exps", type = "string", default = "5,20,35",
         help = "Test n~=2**<these values>."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "b", long = "bits", type = "integer", default = 256L,
         help = "MPFR precision; or if 0, test approx=TRUE and stats against approx=FALSE."),
    list(short = "l", long = "logp", type = "flag",
         help = "Test log.p=TRUE."),
    list(short = "o", long = "omit-stats", type = "flag",
         help = "Skip stats::phyper() calculations."),
    list(short = "s", long = "seed", type = "integer", default = 1L,
         help = "RNG seed.")))
  set.seed(cmd_args$seed)
  phyper_accuracy_test(cmd_args$row_prop,
                       cmd_args$col_prop,
                       cmd_args$z_score,
                       parse_range_string(cmd_args$exps),
                       cmd_args$number,
                       cmd_args$bits,
                       cmd_args$logp,
                       cmd_args$omit_stats)
}

main()
