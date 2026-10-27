#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))
source(file.path(.script_dir, "mpfr_impls.R"))

description <- "Checks the conditional MLE odds ratio and 95% confidence interval returned by
exactr::fisher.test() and (when n < 2^max_stats_exp) stats::fisher.test() on
2x2 tables, against a MPFR-based implementation."

odds_ratio_accuracy_test <- function(ab, ac, z, pow2s, num_trials_per_pow2, bits, max_stats_exp) {
  for (pow2 in pow2s) {
    ns <- sample_pow2_tier(pow2, num_trials_per_pow2)
    relerr_base_est_ssq <- 0
    relerr_base_ci95_low_ssq <- 0
    relerr_base_ci95_high_ssq <- 0
    relerr_stats_est_ssq <- 0
    relerr_stats_ci95_low_ssq <- 0
    relerr_stats_ci95_high_ssq <- 0
    include_stats <- (pow2 < max_stats_exp)
    for (n in ns) {
      tbl <- make_22_table(n, ab, ac, z)
      x <- matrix(tbl, nrow = 2, byrow = TRUE)
      want_est <- cond_odds_ratio(tbl[1], tbl[2], tbl[3], tbl[4], bits)
      want_ci95 <- cond_odds_ratio_ci(tbl[1], tbl[2], tbl[3], tbl[4], bits = bits)
      got_base <- exactr::fisher.test(x)
      relerr <- compute_relerr(got_base$estimate[[1]], want_est)
      relerr_base_est_ssq <- relerr_base_est_ssq + relerr * relerr
      relerr <- compute_relerr(got_base$conf.int[1], want_ci95[1])
      relerr_base_ci95_low_ssq <- relerr_base_ci95_low_ssq + relerr * relerr
      relerr <- compute_relerr(got_base$conf.int[2], want_ci95[2])
      relerr_base_ci95_high_ssq <- relerr_base_ci95_high_ssq + relerr * relerr
      if (include_stats) {
        got_stats <- stats::fisher.test(x)
        relerr <- compute_relerr(got_stats$estimate[[1]], want_est)
        relerr_stats_est_ssq <- relerr_stats_est_ssq + relerr * relerr
        relerr <- compute_relerr(got_stats$conf.int[1], want_ci95[1])
        relerr_stats_ci95_low_ssq <- relerr_stats_ci95_low_ssq + relerr * relerr
        relerr <- compute_relerr(got_stats$conf.int[2], want_ci95[2])
        relerr_stats_ci95_high_ssq <- relerr_stats_ci95_high_ssq + relerr * relerr
      }
    }
    num_trials <- length(ns)
    base_est_rms <- sqrt(relerr_base_est_ssq / num_trials)
    base_ci95_low_rms <- sqrt(relerr_base_ci95_low_ssq / num_trials)
    base_ci95_high_rms <- sqrt(relerr_base_ci95_high_ssq / num_trials)
    print_str <- paste0("n in [2^", pow2, ", 2^", pow2 + 1, "): estErrRMS=", format_g3(base_est_rms),
                        "  ciLowErrRMS=", format_g3(base_ci95_low_rms),
                        "  ciHighErrRMS=", format_g3(base_ci95_high_rms))
    if (include_stats) {
      stats_est_rms <- sqrt(relerr_stats_est_ssq / num_trials)
      stats_ci95_low_rms <- sqrt(relerr_stats_ci95_low_ssq / num_trials)
      stats_ci95_high_rms <- sqrt(relerr_stats_ci95_high_ssq / num_trials)
      print_str <- paste0(print_str, "  statsEstErrRMS=", format_g3(stats_est_rms),
                          "  statsCiLowErrRMS=", format_g3(stats_ci95_low_rms),
                          "  statsCiHighErrRMS=", format_g3(stats_ci95_high_rms))
    }
    cat(print_str, "\n", sep = "")
  }
}

main <- function() {
  cmd_args <- parse_commandline_args(description, list(
    list(short = "z", long = "z-score", type = "numeric", required = TRUE,
         help = "Z-score to test."),
    list(short = "r", long = "row-prop", type = "numeric", default = 0.5,
         help = "Proportion of total in first row."),
    list(short = "c", long = "col-prop", type = "numeric", default = 0.5,
         help = "Proportion of total in first column."),
    list(short = "e", long = "exps", type = "string", default = "5,10,15,35",
         help = "Test n~=2**<these values>."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "b", long = "bits", type = "integer", default = 256L,
         help = "MPFR precision."),
    list(short = "m", long = "max-stats-exp", type = "integer", default = 27L,
         help = "Skip stats::fisher.test() for n >= 2**<this value>."),
    list(short = "s", long = "seed", type = "integer", default = 1L,
         help = "RNG seed.")))
  set.seed(cmd_args$seed)
  odds_ratio_accuracy_test(cmd_args$row_prop,
                           cmd_args$col_prop,
                           cmd_args$z_score,
                           parse_range_string(cmd_args$exps),
                           cmd_args$number,
                           cmd_args$bits,
                           cmd_args$max_stats_exp)
}

main()
