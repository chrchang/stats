#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))

description <- "This collects benchmark results for exactr:::phyper() and stats::phyper()."

phyper_benchmark <- function(ab, ac, z, pow2s, num_trials_per_pow2, logp, omit_stats) {
  warmup <- exactr:::phyper(1, 2, 2, 2)
  if (!omit_stats) {
    warmup <- stats::phyper(1, 2, 2, 2)
  }
  for (pow2 in pow2s) {
    n <- 2^pow2 - 1
    tbl <- make_22_table(n, ab, ac, z)
    a <- tbl[1]
    cur_ab <- tbl[1] + tbl[2]
    cur_ac <- tbl[1] + tbl[3]
    secs_noapprox <- seconds_per_iter(function() exactr:::phyper(a, cur_ac, n - cur_ac, cur_ab, log.p = logp), num_trials_per_pow2)
    secs_approx <- seconds_per_iter(function() exactr:::phyper(a, cur_ac, n - cur_ac, cur_ab, log.p = logp, approx = TRUE), num_trials_per_pow2)
    print_str <- paste0("n=(2^", pow2, ")-1: base=", format_g3(secs_noapprox), "  approx=", format_g3(secs_approx))
    if (!omit_stats) {
      secs_stats <- seconds_per_iter(function() stats::phyper(a, cur_ac, n - cur_ac, cur_ab, log.p = logp), num_trials_per_pow2)
      print_str <- paste0(print_str, "  stats=", format_g3(secs_stats))
    }
    cat(print_str, " sec/iter\n", sep = "")
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
         help = "Test n=2**<these values> - 1."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "l", long = "logp", type = "flag",
         help = "Test log.p=TRUE."),
    list(short = "o", long = "omit-stats", type = "flag",
         help = "Skip stats::phyper() calculations.")))
  phyper_benchmark(cmd_args$row_prop,
                   cmd_args$col_prop,
                   cmd_args$z_score,
                   parse_range_string(cmd_args$exps),
                   cmd_args$number,
                   cmd_args$logp,
                   cmd_args$omit_stats)
}

main()
