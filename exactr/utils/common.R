# Helpers shared by the accuracy and benchmark scripts in this directory.
# These mirror the per-script helpers in Python/utils/; argparse is replaced by
# parse_commandline_args() below since base R has no equivalent.

suppressPackageStartupMessages(library(exactr))

parse_range_string <- function(input_str) {
  result <- c()
  for (part in strsplit(input_str, ",", fixed = TRUE)[[1]]) {
    if (grepl("-", part, fixed = TRUE)) {
      bounds <- as.integer(strsplit(part, "-", fixed = TRUE)[[1]])
      result <- c(result, bounds[1]:bounds[2])
    } else {
      result <- c(result, as.integer(part))
    }
  }
  result
}

flush_if_denormal <- function(x) {
  if (abs(x) < .Machine$double.xmin) {
    return(0)
  }
  x
}

compute_relerr <- function(got, want) {
  if (want == 0) {
    if (got == 0) {
      return(0)
    }
    return(Inf)
  } else if (is.nan(want)) {
    if (is.nan(got)) {
      return(0)
    }
    return(Inf)
  }
  (got - want) / want
}

format_g3 <- function(x) {
  sprintf("%.3g", x)
}

# Equivalent of Python's
#   if num_trials_per_pow2 >= n_limit - min_n:
#     ns = range(min_n, n_limit)
#   else:
#     ns = random.sample(range(min_n, n_limit), k=num_trials_per_pow2)
sample_pow2_tier <- function(pow2, num_trials_per_pow2) {
  min_n <- 2^pow2
  if (num_trials_per_pow2 >= min_n) {
    return(min_n:(2 * min_n - 1))
  }
  min_n - 1 + sample.int(min_n, num_trials_per_pow2)
}

# Returns c(a, b, c, d) for a 2x2 table with total n, first-row total
# ~= ab*n, first-column total ~= ac*n, and upper-left cell at the given z-score.
make_22_table <- function(n, ab, ac, z) {
  cur_ab <- round(ab * n)
  cur_ac <- round(ac * n)
  d_minus_a <- n - cur_ab - cur_ac
  min_a <- max(0, -d_minus_a)
  stdev <- 0
  if (n > 1) {
    variance <- cur_ab * cur_ac * (n - cur_ab) * (n - cur_ac) / (n * n * (n - 1))
    stdev <- sqrt(variance)
  }
  a <- round((cur_ab * cur_ac / n) + z * stdev)
  if (a < min_a) {
    a <- min_a
  } else if (a > min(cur_ab, cur_ac)) {
    a <- min(cur_ab, cur_ac)
  }
  c(a, cur_ab - a, cur_ac - a, d_minus_a + a)
}

seconds_per_iter <- function(f, number) {
  start <- Sys.time()
  for (i in seq_len(number)) {
    f()
  }
  as.numeric(difftime(Sys.time(), start, units = "secs")) / number
}

# Minimal argparse substitute.  spec is a list of options, each of which is a
# list with elements short, long, type ("numeric", "integer", "string", or
# "flag"), default, help, and (optionally) required.  Returns a list keyed by
# long option name, with '-' replaced by '_'.
parse_commandline_args <- function(description, spec) {
  script_name <- basename(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
  print_usage <- function(con) {
    cat("usage: ", script_name, " [options]\n\n", description, "\n\n", sep = "", file = con)
    for (opt in spec) {
      metavar <- if (opt$type == "flag") "" else paste0(" ", toupper(gsub("-", "_", opt$long)))
      default_str <- if (is.null(opt$default) || (opt$type == "flag")) "" else paste0(" (default: ", opt$default, ")")
      required_str <- if (isTRUE(opt$required)) " (required)" else ""
      cat("  -", opt$short, ", --", opt$long, metavar, "\n        ", opt$help, default_str, required_str, "\n", sep = "", file = con)
    }
  }
  fail <- function(...) {
    print_usage(stderr())
    cat("\nerror: ", ..., "\n", sep = "", file = stderr())
    quit(save = "no", status = 2)
  }

  result <- list()
  for (opt in spec) {
    key <- gsub("-", "_", opt$long)
    if (opt$type == "flag") {
      result[[key]] <- FALSE
    } else if (!is.null(opt$default)) {
      result[[key]] <- opt$default
    }
  }
  args <- commandArgs(trailingOnly = TRUE)
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("-h", "--help")) {
      print_usage(stdout())
      quit(save = "no", status = 0)
    }
    value <- NULL
    if (startsWith(arg, "--") && grepl("=", arg, fixed = TRUE)) {
      value <- sub("^[^=]*=", "", arg)
      arg <- sub("=.*$", "", arg)
    }
    opt <- Find(function(o) (arg == paste0("-", o$short)) || (arg == paste0("--", o$long)), spec)
    if (is.null(opt)) {
      fail("unrecognized argument: ", arg)
    }
    key <- gsub("-", "_", opt$long)
    if (opt$type == "flag") {
      result[[key]] <- TRUE
    } else {
      if (is.null(value)) {
        i <- i + 1
        if (i > length(args)) {
          fail("argument ", arg, " expects a value")
        }
        value <- args[i]
      }
      parsed <- switch(opt$type,
                       numeric = suppressWarnings(as.numeric(value)),
                       integer = suppressWarnings(as.integer(value)),
                       string = value)
      if (is.na(parsed)) {
        fail("invalid value for ", arg, ": ", value)
      }
      result[[key]] <- parsed
    }
    i <- i + 1
  }
  for (opt in spec) {
    if (isTRUE(opt$required) && is.null(result[[gsub("-", "_", opt$long)]])) {
      fail("argument --", opt$long, " is required")
    }
  }
  result
}
