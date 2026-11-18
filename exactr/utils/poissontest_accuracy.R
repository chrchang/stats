#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))
source(file.path(.script_dir, "mpfr_impls.R"))

description <- "Checks exactr::poisson.test() and (if log.p=FALSE and r*T < 2^max_stats_exp)
stats::poisson.test()'s one-sample two-sided p-value accuracy against a simple
MPFR-based implementation (slow!).

stats::poisson.test() takes O(r*T) time and memory in the two-sided case, so
it is skipped for larger r*T."

poissontest_accuracy_test <- function(z, pow2s, num_trials_per_pow2, bits, logp, max_stats_exp) {
  for (pow2 in pow2s) {
    relerr_base_ssq <- 0
    relerr_stats_ssq <- 0
    include_stats <- (pow2 < max_stats_exp) && !logp
    for (trial_idx in seq_len(num_trials_per_pow2)) {
      # Non-integer expected count in [2^pow2, 2^(pow2+1)).
      m <- 2^pow2 * (1 + runif(1))
      k <- round(m + z * sqrt(m) + 0.5)
      if (k < 0) {
        k <- 0
      }
      if (logp) {
        got_base <- exactr::poisson.test(k, r = m, log.p = TRUE)$log.p.value
      } else {
        got_base <- exactr::poisson.test(k, r = m)$p.value
      }
      got_base <- flush_if_denormal(got_base)
      want <- flush_if_denormal(poissontest(k, 1, m, bits, logp))
      relerr <- compute_relerr(got_base, want)
      relerr_base_ssq <- relerr_base_ssq + relerr * relerr
      if (include_stats) {
        got_stats <- stats::poisson.test(k, r = m)$p.value
        relerr <- compute_relerr(got_stats, want)
        relerr_stats_ssq <- relerr_stats_ssq + relerr * relerr
      }
    }
    base_rms <- sqrt(relerr_base_ssq / num_trials_per_pow2)
    print_str <- paste0("r*T in [2^", pow2, ", 2^", pow2 + 1, "): errRMS=", format_g3(base_rms))
    if (include_stats) {
      stats_rms <- sqrt(relerr_stats_ssq / num_trials_per_pow2)
      print_str <- paste0(print_str, "  statsErrRMS=", format_g3(stats_rms))
    }
    cat(print_str, "\n", sep = "")
  }
}

main <- function() {
  cmd_args <- parse_commandline_args(description, list(
    # z=0 is not representative.
    list(short = "z", long = "z-score", type = "numeric", required = TRUE,
         help = "Event-count z-score to test."),
    list(short = "e", long = "exps", type = "string", default = "5,20,30",
         help = "Test r*T~=2**<these values>."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "b", long = "bits", type = "integer", default = 256L,
         help = "MPFR precision."),
    list(short = "l", long = "logp", type = "flag",
         help = "Test log.p=TRUE (disables stats)."),
    list(short = "m", long = "max-stats-exp", type = "integer", default = 25L,
         help = "Skip stats::poisson.test() for r*T >= 2**<this value>."),
    list(short = "s", long = "seed", type = "integer", default = 1L,
         help = "RNG seed.")))
  set.seed(cmd_args$seed)
  poissontest_accuracy_test(cmd_args$z_score,
                            parse_range_string(cmd_args$exps),
                            cmd_args$number,
                            cmd_args$bits,
                            cmd_args$logp,
                            cmd_args$max_stats_exp)
}

main()
