###############################################################################
# 99_verify_all.R  --  single, self-contained verification script (base R only)
#
# Run it like the other scripts (Rscript 99_verify_all.R, or RStudio "Source"): data are
# read from DATA_DIR and verify_closed_forms.csv is written to OUT_DIR (paths.R).
# No other script, cache or package is needed. It is not part of run_local.bat.
#
#   Part 1  data files : row/column counts of the five CSVs
#   Part 2  closed forms of Theorem 5 and Supplementary Note D (Appendix D)
#           against the recursion of Corollary 1
#           (E(R), Var(R), Cov(tau,R), Var(tau); n = 256, 512; 21 finite d)
#   Part 3  implementation checks (Supplementary Note 7)
#           (a) absorption (Lemma 5), (b) node-counting identity (Lemma 6),
#           (c) kernel rows of Theorem 4 sum to one
#
# Expected: Part 2 max difference below 1e-13; Part 3 all violations 0 and
# kernel deviations at floating-point level (~1e-16). Run time: a few minutes.
###############################################################################
rm(list = ls())
t_start <- Sys.time()

script_dir <- function() {
  a <- commandArgs(FALSE); k <- grep("^--file=", a)
  if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]]))))
  for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }
  getwd()
}
source(file.path(Sys.getenv("CODE_DIR", script_dir()), "paths.R"))   # DATA_DIR, OUT_DIR
cat("\n")

Hn  <- function(n) if (n < 1) 0 else sum(1 / seq_len(n))
Hn2 <- function(n) if (n < 1) 0 else sum(1 / seq_len(n)^2)

###############################################################################
## Part 1  data files
###############################################################################
cat("== Part 1: data files ==\n")
expect <- c(quakes = 5, magic = 10, wine_white = 11, letter = 16, supercond = 81)
for (nm in names(expect)) {
  f <- file.path(DATA_DIR, paste0(nm, ".csv"))
  if (!file.exists(f)) { cat(sprintf("  %-15s MISSING\n", basename(f))); next }
  x <- read.csv(f)
  cat(sprintf("  %-15s %6d rows x %3d cols  (expected %2d cols)  %s\n", basename(f), nrow(x), ncol(x),
              expect[[nm]], if (ncol(x) == expect[[nm]]) "OK" else "CHECK"))
}

###############################################################################
## Part 2  closed forms vs recursion
###############################################################################
## ---- recursion of Corollary 1 (same as 01_exact_simul_figures.R) -------------
bern_sum <- function(ps) { p <- 1; for (q in ps) p <- c(p * (1 - q), 0) + c(0, p * q); p }
ext_dists <- function(n) lapply(seq_len(n), function(m) if (m == 1) 1 else bern_sum(1 / seq_len(m - 1)))

kernels <- function(n, d) {                     # Theorem 4: rows m, columns j
  Kabs <- matrix(0, n, n); KG <- matrix(0, n, n)
  for (m in 3:n) {
    j <- seq_len(m - 1)
    pm <- (2 / m)^d; am <- (2 / m)^(d - 1); den <- m * (m - 1) * (1 - pm)
    Kabs[m, 1] <- 2 * (1 - am) / den
    Kabs[m, 2:(m - 1)] <- 2 * (2 * (2 / j[-1])^(d - 1) - am) / den
    if (m > 3) KG[m, 3:(m - 1)] <- 2 * (j[3:(m - 1)] - 2 * (2 / j[3:(m - 1)])^(d - 1)) / den
  }
  list(Kabs = Kabs, KG = KG)
}
joint_tau_M <- function(n, d) {                 # J[t+1, m] = P(tau = t, M_tau = m)
  K <- kernels(n, d); pn <- (2 / n)^d
  J <- matrix(0, n + 1, n); f <- numeric(n); J[1, n] <- pn; f[n] <- 1 - pn; h <- 0
  while (sum(f) > 1e-16 && h < n) {
    J[h + 2, ] <- J[h + 2, ] + as.vector(f %*% K$Kabs)
    f <- as.vector(f %*% K$KG); h <- h + 1
  }
  J[seq_len(max(which(rowSums(J) > 0))), , drop = FALSE]
}
recursion_moments <- function(n, d, ext) {
  J <- joint_tau_M(n, d); P <- matrix(0, nrow(J), n)
  for (m in seq_len(n)) { w <- J[, m]; if (!any(w > 0)) next
    e <- ext[[m]]; P[, seq_along(e)] <- P[, seq_along(e)] + outer(w, e) }
  t <- 0:(nrow(P) - 1); r <- 0:(ncol(P) - 1); pt <- rowSums(P); pr <- colSums(P)
  Et <- sum(t * pt); ER <- sum(r * pr)
  c(E_R = ER, Var_R = sum(r^2 * pr) - ER^2, Cov = sum(outer(t, r) * P) - Et * ER,
    Var_tau = sum(t^2 * pt) - Et^2, E_H = Et + ER)
}

## ---- closed forms (Theorem 5; Supplementary Note D) --------------------------
## sums over every internal node size m = 2..n, the root (m = n) included:
##   c_n(m) = 2n/{m(m+1)} for m < n, c_n(n) = 1          (Lemma 7)
##   s(m)   = 2^d / m^(d-1)                                (extreme observations in a node of size m)
##   h(m)   = H_{m-1}                                      (Theorem 1)
##   D_n(m) = 2(H_n - H_{m+1}) + 1 for m < n, D_n(n) = 0   (Appendix C.2; the root is at depth 0)
closed_forms <- function(n, d) {
  m  <- 2:n
  cn <- ifelse(m == n, 1, 2 * n / (m * (m + 1)))
  s  <- 2^d / m^(d - 1)
  h  <- sapply(m, function(k) Hn(k - 1))
  Dn <- sapply(m, function(k) if (k == n) 0 else 2 * (Hn(n) - Hn(k + 1)) + 1)
  E_R   <- sum(cn * s) / n                         # Theorem 5
  E_R2  <- sum(cn * s * (2 * h - 1)) / n           # Note D.1
  Var_R <- E_R2 - E_R^2
  E_HR  <- sum(cn * s * (Dn + h)) / n              # Note D.3
  E_H   <- 2 * Hn(n) - 2; Var_H <- 2 * Hn(n) + 2 - 4 * Hn2(n)   # Theorem 3
  Cov   <- E_HR - E_H * E_R - Var_R                # Note D.2
  c(E_R = E_R, Var_R = Var_R, Cov = Cov, Var_tau = Var_H - Var_R - 2 * Cov)   # Note D.3
}

cat("\n== Part 2: closed forms vs recursion ==\n")
d_set <- sort(unique(c(1:10, round(10^seq(1, 2, length.out = 12)))))
rows <- list()
for (n in c(256L, 512L)) {
  ext <- ext_dists(n)
  for (d in d_set) {
    rc <- recursion_moments(n, d, ext); cf <- closed_forms(n, d)
    q <- names(cf)
    rows[[length(rows) + 1]] <- data.frame(n = n, d = d, quantity = q, closed_form = unname(cf),
                                           recursion = unname(rc[q]), abs_diff = abs(unname(cf) - unname(rc[q])))
    if (d == 5) cat(sprintf("  n = %d: |E(H) - (2H_n - 2)| from the recursion = %.2e\n", n, abs(rc[["E_H"]] - (2 * Hn(n) - 2))))
  }
}
CF <- do.call(rbind, rows)
write.csv(CF, file.path(OUT_DIR, "verify_closed_forms.csv"), row.names = FALSE)
cat(sprintf("  compared %d values (n = 256, 512; %d finite d; 4 quantities)\n", nrow(CF), length(d_set)))
print(aggregate(abs_diff ~ quantity, CF, max))
cat(sprintf("  MAXIMUM absolute difference: %.3e\n", max(CF$abs_diff)))

###############################################################################
## Part 3  implementation checks
###############################################################################
cat("\n== Part 3: implementation checks ==\n")
is_extreme <- function(S, k) {                  # minimum or maximum in every column of S
  row <- matrix(S[k, ], nrow(S), ncol(S), byrow = TRUE)
  all(colSums(S < row) == 0L | colSums(S > row) == 0L)
}
one_tree_checked <- function(Z, d) {            # splitting rule of Methods, instrumented
  n <- nrow(Z); H <- integer(n); TAU <- rep(NA_integer_, n)
  ext_now <- logical(n); s_total <- 0L; viol <- 0L
  stack <- list(list(idx = seq_len(n), h = 0L)); top <- 1L
  while (top > 0L) {
    nd <- stack[[top]]; stack[[top]] <- NULL; top <- top - 1L
    idx <- nd$idx; h <- nd$h; m <- length(idx)
    if (m == 1L) { H[idx] <- h; if (is.na(TAU[idx])) TAU[idx] <- h; next }
    S <- Z[idx, , drop = FALSE]
    now <- vapply(seq_len(m), function(k) is_extreme(S, k), logical(1))
    s_total <- s_total + sum(now)                                   # s(A)
    viol <- viol + sum(ext_now[idx] & !now)                          # (a)
    ext_now[idx] <- now
    newly <- idx[now & is.na(TAU[idx])]; TAU[newly] <- h
    repeat { q <- sample.int(d, 1L); v <- S[, q]; rg <- range(v); if (rg[2] > rg[1]) break }
    cut <- runif(1L, rg[1], rg[2])
    top <- top + 1L; stack[[top]] <- list(idx = idx[v < cut],  h = h + 1L)
    top <- top + 1L; stack[[top]] <- list(idx = idx[v >= cut], h = h + 1L)
  }
  c(sum_R = sum(H - TAU), sum_s = s_total, absorption_viol = viol)
}
n_chk <- 64L; n_trees <- 200L
for (d in c(1L, 2L, 3L, 5L)) {
  set.seed(20260929L + d)
  res <- t(replicate(n_trees, one_tree_checked(matrix(runif(n_chk * d), n_chk, d), d)))
  cat(sprintf("  (a)(b) d = %d: node-counting identity violated in %d/%d trees; absorption violations: %d\n",
              d, sum(res[, "sum_R"] != res[, "sum_s"]), n_trees, sum(res[, "absorption_viol"])))
}
for (n in c(32L, 128L)) for (d in c(1L, 2L, 5L, 10L)) {
  K <- kernels(n, d); dev <- max(abs(rowSums(K$Kabs + K$KG)[3:n] - 1))
  cat(sprintf("  (c) n = %3d, d = %2d: max |kernel row sum - 1| = %.2e\n", n, d, dev))
}

cat(sprintf("\nwritten: %s\nelapsed: %.1f min\n", file.path(OUT_DIR, "verify_closed_forms.csv"),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
