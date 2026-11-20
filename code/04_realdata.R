###############################################################################
# Runs the real-data pipeline (Results: real data) for all five data sets and both
# marginals ("rank" and "raw"), writing one <name>_<marginal>_components.csv
# per combination to /results. These files are the input to 05_table2_extdata_table2.R.
#
# Data required in /data (explanatory variables only, no target column):
#   quakes.csv, magic.csv, wine_white.csv, letter.csv, supercond.csv
###############################################################################
rm(list = ls())
suppressWarnings(suppressMessages(library(parallel)))
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "04_realdata_core.R"))

CFG <- list(
  datasets  = c("quakes", "magic", "wine_white", "letter", "supercond"),
  marginals = c("rank", "raw"),
  psi       = 256L,
  n_sub     = if (FULL_RUN) 20L else 8L,     # FULL=1: 20 subsamples (paper); default: 8 (quick check)
  n_trees   = if (FULL_RUN) 200L else 100L,
  jitter_sd = 1e-6,
  seed      = 20260922L,
  cores     = max(1L, parallel::detectCores() - 1L)
)

for (ds in CFG$datasets) for (mg in CFG$marginals) {
  cat(sprintf("\n== %s / %s ==\n", ds, mg))
  out <- run_one_dataset(ds, mg, DATA_DIR, CFG$psi, CFG$n_sub, CFG$n_trees, CFG$jitter_sd, CFG$seed, CFG$cores)
  outf <- file.path(OUT_DIR, sprintf("realdata_%s_%s.csv", ds, mg))
  write.csv(rbind(cbind(kind = "moment", out$moment), cbind(kind = "D", out$D)), outf, row.names = FALSE)
  cat(sprintf("N = %d, d = %d, psi = %d, %d subsamples x %d trees -> %s\n",
              out$N, out$d, CFG$psi, CFG$n_sub, CFG$n_trees, outf))
}
cat("\nAll components written to", OUT_DIR, "\n")
