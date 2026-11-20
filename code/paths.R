## Paths shared by all scripts.
## On Code Ocean: /data is read-only input, /results is the output folder.
## Elsewhere: data are read from LOCAL_DATA and results written to LOCAL_OUT when
## these environment variables are set (run_local.bat sets them); otherwise the
## repository layout is assumed, i.e. ../data and ../results relative to this file.
on_capsule <- dir.exists("/data") && dir.exists("/results") && dir.exists("/code")

.paths_dir <- local({                     # folder that holds paths.R (works for Rscript and RStudio "Source")
  f <- NULL
  for (i in rev(seq_len(sys.nframe()))) { o <- sys.frame(i)$ofile; if (!is.null(o)) { f <- o; break } }
  if (is.null(f)) Sys.getenv("CODE_DIR", getwd()) else dirname(normalizePath(f))
})
.repo_root <- normalizePath(file.path(.paths_dir, ".."), mustWork = FALSE)

DATA_DIR <- if (on_capsule) "/data"    else Sys.getenv("LOCAL_DATA", file.path(.repo_root, "data"))
OUT_DIR  <- if (on_capsule) "/results" else Sys.getenv("LOCAL_OUT",  file.path(.repo_root, "results"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
FULL_RUN <- identical(Sys.getenv("FULL", "0"), "1")   # FULL=1: recompute everything at the paper's settings
cat("DATA_DIR =", DATA_DIR, "\nOUT_DIR  =", OUT_DIR, "\nFULL_RUN =", FULL_RUN, "\n")
if (!dir.exists(DATA_DIR)) stop("data folder not found: ", DATA_DIR, " (set LOCAL_DATA or keep the repository layout)")
