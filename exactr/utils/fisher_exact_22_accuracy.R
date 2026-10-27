#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))
source(file.path(.script_dir, "mpfr_impls.R"))

description <- "Checks exactr::fisher.test() and (when log.p=FALSE and n < 2^max_stats_exp)
stats::fisher.test()'s two-sided 2x2 p-value accuracy against a MPFR-based
implementation (slow!).

stats::fisher.test() requires 32-bit integer table entries, and allocates O(n)
memory, so it is skipped for larger n."

fisher_exact_22_accuracy_test <- function(ab, ac, z, pow2s, num_trials_per_pow2, bits, logp, max_stats_exp) {
  for (pow2 in pow2s) {
    ns <- sample_pow2_tier(pow2, num_trials_per_pow2)
    relerr_base_ssq <- 0
    relerr_stats_ssq <- 0
    include_stats <- (pow2 < max_stats_exp) && !logp
    for (n in ns) {
      tbl <- make_22_table(n, ab, ac, z)
      x <- matrix(tbl, nrow = 2, byrow = TRUE)
      want <- flush_if_denormal(fisher_exact_22(tbl[1], tbl[2], tbl[3], tbl[4], bits, logp))
      if (logp) {
        got_base <- exactr::fisher.test(x, conf.int = FALSE, log.p = TRUE)$log.p.value
      } else {
        got_base <- exactr::fisher.test(x, conf.int = FALSE)$p.value
      }
      got_base <- flush_if_denormal(got_base)
      relerr <- compute_relerr(got_base, want)
      relerr_base_ssq <- relerr_base_ssq + relerr * relerr
      if (include_stats) {
        got_stats <- stats::fisher.test(x, conf.int = FALSE)$p.value
        relerr <- compute_relerr(got_stats, want)
        relerr_stats_ssq <- relerr_stats_ssq + relerr * relerr
      }
    }
    num_trials <- length(ns)
    base_rms <- sqrt(relerr_base_ssq / num_trials)
    print_str <- paste0("n in [2^", pow2, ", 2^", pow2 + 1, "): errRMS=", format_g3(base_rms))
    if (include_stats) {
      stats_rms <- sqrt(relerr_stats_ssq / num_trials)
      print_str <- paste0(print_str, "  statsErrRMS=", format_g3(stats_rms))
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
    list(short = "e", long = "exps", type = "string", default = "5,20,35",
         help = "Test n~=2**<these values>."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "b", long = "bits", type = "integer", default = 256L,
         help = "MPFR precision."),
    list(short = "l", long = "logp", type = "flag",
         help = "Test log.p=TRUE (disables stats)."),
    list(short = "m", long = "max-stats-exp", type = "integer", default = 27L,
         help = "Skip stats::fisher.test() for n >= 2**<this value>."),
    list(short = "s", long = "seed", type = "integer", default = 1L,
         help = "RNG seed.")))
  set.seed(cmd_args$seed)
  fisher_exact_22_accuracy_test(cmd_args$row_prop,
                                cmd_args$col_prop,
                                cmd_args$z_score,
                                parse_range_string(cmd_args$exps),
                                cmd_args$number,
                                cmd_args$bits,
                                cmd_args$logp,
                                cmd_args$max_stats_exp)
}

main()
