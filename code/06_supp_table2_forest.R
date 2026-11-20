###############################################################################
# Forest-level decomposition of depth variance   (Supplementary Note 8, Supplementary Table 2)
#
# For a fixed sample X, each observation i has a forest mean depth
#     mu_i = E_T(H_i | X)      (limit of the forest average as T -> Inf)
# and a tree-level variance
#     s2_i = Var_T(H_i | X).
# The single-tree variance splits as
#     Var_{i,T}(H | X) = Var_i(mu_i) + E_i(s2_i)
#            (between observations) + (between trees).
# The same decomposition is computed for tau (entry depth) and R = H - tau.
#
# Output (next to this script):
#   IFforest_decomp.csv   one row per (n, d): the components, averaged over samples,
#                         with standard errors across samples, and the
#                         between-sample variance of the sample mean
###############################################################################

rm(list = ls())
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))
suppressWarnings(suppressMessages(library(parallel)))

QUICK <- !FULL_RUN           # quick (8 x 200) unless FULL=1 is set (20 x 500)
CFG <- list(
  n       = 256L,
  dims    = c(1L, 2L, 3L, 5L, 10L, 20L, 50L, 100L),
  n_trees = if (QUICK) 200L else 500L,   # T per sample
  n_samp  = if (QUICK) 8L   else 20L,    # independent samples per (n, d)
  seed    = 20260921L,
  cores   = max(1L, parallel::detectCores() - 1L)
)

message("output folder: ", OUT_DIR)

# ---------- one tree: returns H, tau for every observation (the splitting rule of Methods)
one_tree <- function(X) {
  n <- nrow(X); d <- ncol(X)
  H <- integer(n); TAU <- rep(NA_integer_, n)
  # stack of nodes: list(idx, depth)
  stack <- vector("list", 2L * n); top <- 0L
  top <- top + 1L; stack[[top]] <- list(idx = seq_len(n), h = 0L)
  while (top > 0L) {
    nd <- stack[[top]]; top <- top - 1L
    idx <- nd$idx; h <- nd$h; m <- length(idx)
    if (m == 1L) { H[idx] <- h; if (is.na(TAU[idx])) TAU[idx] <- h; next }
    S <- X[idx, , drop = FALSE]
    # extreme state: min or max in every coordinate. A point extreme everywhere
    # must be the min or the max of coordinate 1, so at most two candidates.
    v1 <- S[, 1L]
    cand <- idx[v1 == min(v1) | v1 == max(v1)]
    cand <- cand[is.na(TAU[cand])]
    for (cc in cand) {
      row <- rep(X[cc, ], each = m)              # column-major layout of an m x d matrix
      below <- colSums(S < row); above <- colSums(S > row)
      if (all(below == 0L | above == 0L)) TAU[cc] <- h
    }
    q <- sample.int(d, 1L)
    v <- S[, q]; rg <- range(v)
    cut <- runif(1L, rg[1L], rg[2L])
    L <- idx[v < cut]; R <- idx[v >= cut]
    if (length(L) == 0L || length(R) == 0L) { top <- top + 1L; stack[[top]] <- nd; next }
    top <- top + 1L; stack[[top]] <- list(idx = L, h = h + 1L)
    top <- top + 1L; stack[[top]] <- list(idx = R, h = h + 1L)
  }
  list(H = H, tau = TAU)
}

# ---------- one sample: T trees -> per-observation mean and variance over trees
one_sample <- function(n, d, n_trees, seed) {
  set.seed(seed)
  X <- matrix(runif(n * d), n, d)
  Hm <- matrix(0L, n_trees, n); Tm <- matrix(0L, n_trees, n)
  for (t in seq_len(n_trees)) { r <- one_tree(X); Hm[t, ] <- r$H; Tm[t, ] <- r$tau }
  Rm <- Hm - Tm
  decomp <- function(M) {
    mu <- colMeans(M)                          # mu_i
    s2 <- colMeans(M^2) - mu^2                 # Var_T(. | i)   (population)
    tot <- mean(M^2) - mean(M)^2               # Var_{i,T}(. | X)
    c(total = tot, between_obs = mean(mu^2) - mean(mu)^2, between_tree = mean(s2), mean = mean(M))
  }
  c(H = decomp(Hm), tau = decomp(Tm), R = decomp(Rm))
}

# ---------- run 
rows <- list()
for (d in CFG$dims) {
  message(sprintf("n = %d, d = %d : %d samples x %d trees", CFG$n, d, CFG$n_samp, CFG$n_trees))
  seeds <- CFG$seed + 1000L * d + seq_len(CFG$n_samp)
  res <- if (.Platform$OS.type == "windows") {
    cl <- makeCluster(CFG$cores)
    clusterExport(cl, c("one_sample", "one_tree", "CFG"))
    r <- parLapply(cl, seeds, function(s, n, d, n_trees) one_sample(n, d, n_trees, s),
                   n = CFG$n, d = d, n_trees = CFG$n_trees)
    stopCluster(cl); r
  } else mclapply(seeds, function(s) one_sample(CFG$n, d, CFG$n_trees, s), mc.cores = CFG$cores)
  M <- do.call(rbind, res)                       # n_samp x 12
  mean_over_samples <- colMeans(M)
  se_over_samples   <- apply(M, 2L, sd) / sqrt(nrow(M))
  # between-sample variance of the sample mean (the term that Var(H) adds on top of E_X Var(H|X))
  bsv <- c(H = var(M[, "H.mean"]), tau = var(M[, "tau.mean"]), R = var(M[, "R.mean"]))
  row <- data.frame(n = CFG$n, d = d, t(mean_over_samples), t(setNames(se_over_samples, paste0(names(se_over_samples), "_se"))),
                    between_sample_H = bsv[["H"]], between_sample_tau = bsv[["tau"]], between_sample_R = bsv[["R"]])
  rows[[length(rows) + 1L]] <- row
  write.csv(do.call(rbind, rows), file.path(OUT_DIR, "SuppTable2_forest_variance.csv"), row.names = FALSE)
}
RES <- do.call(rbind, rows)

Hn <- function(n) sum(1 / seq_len(n)); Hn2 <- function(n) sum(1 / seq_len(n)^2)
theory_VarH <- 2 * Hn(CFG$n) + 2 - 4 * Hn2(CFG$n)
cat(sprintf("\nn = %d, theory Var(H) = %.4f (dimension-free)\n", CFG$n, theory_VarH))
cat("  d   Var(H|X)  between_obs  between_tree |  Var(tau|X) between_obs between_tree |  Var(R|X) between_obs between_tree | share_obs(H)\n")
for (i in seq_len(nrow(RES))) with(RES[i, ], cat(sprintf(
  "%3d   %7.3f   %8.3f     %8.3f   |   %7.3f   %8.3f    %8.3f   |  %7.3f   %8.3f    %8.3f   |   %5.3f\n",
  d, H.total, H.between_obs, H.between_tree, tau.total, tau.between_obs, tau.between_tree,
  R.total, R.between_obs, R.between_tree, H.between_obs / H.total)))
cat("\nwritten:", file.path(OUT_DIR, "SuppTable2_forest_variance.csv"), "\n")
