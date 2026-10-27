#!/usr/bin/env Rscript
.script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
source(file.path(.script_dir, "common.R"))

description <- "This collects 2x2 benchmark results for exactr::fisher.test() and
stats::fisher.test() with conf.int=TRUE (the default), i.e. p-value +
conditional MLE odds ratio + 95% confidence interval.  stats::fisher.test() is
skipped for n >= 2^max_stats_exp."

odds_ratio_benchmark <- function(ab, ac, z, pow2s, num_trials_per_pow2, max_stats_exp) {
  warmup <- exactr::fisher.test(matrix(1, 2, 2))
  warmup <- stats::fisher.test(matrix(1, 2, 2))
  for (pow2 in pow2s) {
    n <- 2^pow2 - 1
    x <- matrix(make_22_table(n, ab, ac, z), nrow = 2, byrow = TRUE)
    secs_base <- seconds_per_iter(function() exactr::fisher.test(x), num_trials_per_pow2)
    print_str <- paste0("n=(2^", pow2, ")-1: base=", format_g3(secs_base))
    if (pow2 < max_stats_exp) {
      secs_stats <- seconds_per_iter(function() stats::fisher.test(x), num_trials_per_pow2)
      print_str <- paste0(print_str, "  stats=", format_g3(secs_stats))
    }
    cat(print_str, " sec/iter\n", sep = "")
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
         help = "Test n=2**<these values> - 1."),
    list(short = "n", long = "number", type = "integer", default = 10L,
         help = "Number of trials per power-of-2 tier."),
    list(short = "m", long = "max-stats-exp", type = "integer", default = 27L,
         help = "Skip stats::fisher.test() for n >= 2**<this value>.")))
  odds_ratio_benchmark(cmd_args$row_prop,
                       cmd_args$col_prop,
                       cmd_args$z_score,
                       parse_range_string(cmd_args$exps),
                       cmd_args$number,
                       cmd_args$max_stats_exp)
}

main()
