###############################################################################
# Supplementary Tables 4 and 5 (Supplementary Note 9)
#
# Supp. Table 4  number of variables cut at once.  At every node choose r distinct
#        variables, draw one cut point per chosen variable uniformly between the
#        node's min and max of that variable, and split the node into the
#        non-empty cells of the resulting 2^r-way partition.  r = 1 is the model
#        of Methods.  r = d cuts all variables (quadtree / octree).
#        Output: E(H), Var(H) of a uniformly chosen observation, by (r, d).
#
# Supp. Table 5  node occupancy under binary vs 4-ary splits on one variable.  For every
#        internal node of size m, count the observations that are extreme (min or
#        max) in every coordinate; compare the mean count with 2^d / m^(d-1).
#
# Output (next to this script): IFS3_split_rules.csv, IFS3_occupancy.csv
# base R + parallel.  Windows uses a socket cluster.
###############################################################################

rm(list = ls())
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))
suppressWarnings(suppressMessages(library(parallel)))

QUICK <- !FULL_RUN
CFG <- list(
  n       = 256L,
  dims    = c(2L, 3L, 5L, 10L),
  r_fixed = c(1L, 2L, 3L),
  n_trees = if (QUICK) 200L else 500L,
  n_rep   = if (QUICK) 8L   else 20L,
  occ_d   = 3L,
  occ_m   = c(2L, 4L, 8L, 16L),
  seed    = 20260921L,
  cores   = max(1L, parallel::detectCores() - 1L)
)

## ---------------------------------------------------------------------------
## Supp. Table 4: tree with r variables cut per split; returns depth of every observation
## ---------------------------------------------------------------------------
tree_r <- function(X, r) {
  n <- nrow(X); d <- ncol(X); H <- integer(n)
  stack <- vector("list", 4L * n); top <- 1L; stack[[1L]] <- list(idx = seq_len(n), h = 0L)
  while (top > 0L) {
    nd <- stack[[top]]; top <- top - 1L; idx <- nd$idx; h <- nd$h; m <- length(idx)
    if (m == 1L) { H[idx] <- h; next }
    qs <- if (r >= d) seq_len(d) else sample.int(d, r)
    S <- X[idx, qs, drop = FALSE]
    key <- integer(m)
    for (k in seq_along(qs)) { v <- S[, k]; rg <- range(v); cut <- runif(1L, rg[1L], rg[2L])
      key <- key + (v >= cut) * 2L^(k - 1L) }
    cells <- split(idx, key)
    if (length(cells) == 1L) { top <- top + 1L; stack[[top]] <- nd; next }   # degenerate, redo
    for (cl in cells) { top <- top + 1L; stack[[top]] <- list(idx = cl, h = h + 1L) }
  }
  H
}
rep_split <- function(n, d, r, n_trees, seed) {
  set.seed(seed); X <- matrix(runif(n * d), n, d)
  Hs <- unlist(lapply(seq_len(n_trees), function(t) tree_r(X, r)))
  c(mean_H = mean(Hs), var_H = mean(Hs^2) - mean(Hs)^2)
}

## ---------------------------------------------------------------------------
## Supp. Table 5: occupancy s(A) by node size under binary and 4-ary splits (one variable)
## ---------------------------------------------------------------------------
occupancy_tree <- function(X, arity, m_list) {
  n <- nrow(X); d <- ncol(X)
  cnt <- setNames(numeric(length(m_list)), m_list); tot <- cnt
  stack <- vector("list", 4L * n); top <- 1L; stack[[1L]] <- seq_len(n)
  while (top > 0L) {
    idx <- stack[[top]]; top <- top - 1L; m <- length(idx)
    if (m == 1L) next
    S <- X[idx, , drop = FALSE]
    if (m %in% m_list) {
      # extreme in every coordinate: candidates are the min and max of coordinate 1
      v1 <- S[, 1L]; cand <- which(v1 == min(v1) | v1 == max(v1)); s <- 0L
      for (cc in cand) { row <- rep(S[cc, ], each = m)
        if (all(colSums(S < row) == 0L | colSums(S > row) == 0L)) s <- s + 1L }
      key <- as.character(m); cnt[key] <- cnt[key] + s; tot[key] <- tot[key] + 1
    }
    q <- sample.int(d, 1L); v <- S[, q]; rg <- range(v)
    cuts <- sort(runif(arity - 1L, rg[1L], rg[2L]))
    key <- findInterval(v, cuts)
    cells <- split(idx, key)
    if (length(cells) == 1L) { top <- top + 1L; stack[[top]] <- idx; next }
    for (cl in cells) { top <- top + 1L; stack[[top]] <- cl }
  }
  rbind(cnt, tot)
}
rep_occ <- function(n, d, arity, m_list, n_trees, seed) {
  set.seed(seed); X <- matrix(runif(n * d), n, d)
  acc <- Reduce(`+`, lapply(seq_len(n_trees), function(t) occupancy_tree(X, arity, m_list)))
  acc[1, ] / pmax(acc[2, ], 1)                                 # mean s(A) by size
}

## ---------------------------------------------------------------------------
## run helpers (parallel over replications)
## ---------------------------------------------------------------------------
run_par <- function(fun, seeds, ...) {
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(CFG$cores); on.exit(stopCluster(cl))
    clusterExport(cl, c("tree_r", "rep_split", "occupancy_tree", "rep_occ", "CFG"))
    parLapply(cl, seeds, fun, ...)
  } else mclapply(seeds, fun, ..., mc.cores = CFG$cores)
}

## Supp. Table 4 ---------------------------------------------------------------------
rows <- list()
for (d in CFG$dims) for (r in unique(c(CFG$r_fixed[CFG$r_fixed <= d], d))) {
  message(sprintf("Supp. Table 4  d = %d, r = %d%s", d, r, if (r == d) " (all variables)" else ""))
  seeds <- CFG$seed + 100L * d + 10L * r + seq_len(CFG$n_rep)
  res <- run_par(function(s, n, d, r, nt) rep_split(n, d, r, nt, s), seeds,
                 n = CFG$n, d = d, r = r, nt = CFG$n_trees)
  M <- do.call(rbind, res)
  rows[[length(rows) + 1L]] <- data.frame(n = CFG$n, d = d, r = r, all_variables = (r == d),
    mean_H = mean(M[, "mean_H"]), mean_H_se = sd(M[, "mean_H"]) / sqrt(nrow(M)),
    var_H = mean(M[, "var_H"]), var_H_se = sd(M[, "var_H"]) / sqrt(nrow(M)))
  write.csv(do.call(rbind, rows), file.path(OUT_DIR, "SuppTable4_split_rules.csv"), row.names = FALSE)
}
SPLIT <- do.call(rbind, rows)
Hn <- sum(1 / seq_len(CFG$n))
cat(sprintf("\nSupp. Table 4  E(H) by r (rows) and d (cols); r = 1 theory = 2H_n - 2 = %.4f\n", 2 * Hn - 2))
print(reshape(SPLIT[, c("d", "r", "mean_H")], idvar = "r", timevar = "d", direction = "wide"), row.names = FALSE)

## Supp. Table 5 ---------------------------------------------------------------------
occ <- list()
for (arity in c(2L, 4L)) {
  message(sprintf("Supp. Table 5  arity = %d", arity))
  seeds <- CFG$seed + 1000L * arity + seq_len(CFG$n_rep)
  res <- run_par(function(s, n, d, ar, ml, nt) rep_occ(n, d, ar, ml, nt, s), seeds,
                 n = CFG$n, d = CFG$occ_d, ar = arity, ml = CFG$occ_m, nt = CFG$n_trees)
  M <- do.call(rbind, res)
  occ[[as.character(arity)]] <- colMeans(M)
}
OCC <- data.frame(m = CFG$occ_m, binary = occ[["2"]], four_ary = occ[["4"]],
                  theory = 2^CFG$occ_d / CFG$occ_m^(CFG$occ_d - 1))
write.csv(OCC, file.path(OUT_DIR, "SuppTable5_occupancy.csv"), row.names = FALSE)
cat(sprintf("\nSupp. Table 5  mean number of extreme observations in a node of size m (d = %d)\n", CFG$occ_d))
print(OCC, row.names = FALSE, digits = 4)
cat("\nwritten:", file.path(OUT_DIR, "SuppTable4_split_rules.csv"), "and SuppTable5_occupancy.csv\n")
