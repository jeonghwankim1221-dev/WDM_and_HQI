###############################################################################
# 10_supp_tables7_9_maximal_state.R -- Supplementary Note 12 (maximal state): Supplementary Tables 7-9
#
#   Table 7 : mean entry depth into the sigma-maximal state, E(tau'_sigma), exact from
#             Supplementary Equation (41) with the recursion of Bentley et al.,
#             A_d(m) = A_d(m-1) + A_{d-1}(m)/m, A_1 = 1;
#             last row: lim_n E(R'_sigma) = 2(d-1)zeta(d) - 1.
#   Table 8 : division of the depth at entry into maximality in some orientation (n = 256).
#             E(tau'_U) by Corollary S4 with B_d(m) = m for m <= 2^d and, for m > 2^d,
#             by Monte Carlo from max(40, floor(6000/m)) samples of size m; E(R) exact (Theorem 5).
#   Table 9 : mean and variance decomposition at entry into the sigma-maximal state (n = 256),
#             tree simulation, N9 samples x 100 trees (environment variable N9, default 50).
#
#   Basic model: n iid Uniform observations on [0,1]^d; at each split one variable is chosen
#   uniformly and the cut point is uniform between its minimum and maximum in the node.
#   Orientation sigma = (+,...,+): X_i is sigma-maximal in a node if no other observation of the
#   node exceeds it in every variable.
#   Outputs: OUT_DIR/maximal_table7.csv, maximal_table8.csv, maximal_table9.csv
#   Quick check of the script: set N9=2 (Table 9 with 2 samples) before the full run.
###############################################################################
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))
library(parallel)

N9    <- as.integer(Sys.getenv("N9", "50"))
SEED  <- 20261003L
CORES <- max(1L, detectCores() - 1L)

harmonic <- function(n, p = 1) sum(1 / seq_len(n)^p)
c_n      <- function(n, m) ifelse(m == n, 1, 2 * n / (m * (m + 1)))          # Lemma 7
zeta     <- function(s, K = 200000) sum(1 / seq_len(K)^s)

A_table <- function(dmax, M) {                                              # A[d, m], m = 1..M
  A <- matrix(0, dmax, M)
  A[1, ] <- 1
  for (d in 2:dmax) for (m in seq_len(M))
    A[d, m] <- (if (m == 1) 0 else A[d, m - 1]) + A[d - 1, m] / m
  A
}
E_Rprime <- function(n, d, A) { m <- 2:n; sum(c_n(n, m) * A[d, m]) / n }   # Proposition S6

## ------------------------------------------------------------------ Table 7 (exact)
D7 <- c(1, 2, 3, 5, 8, 10, 15, 20, 30); N7 <- c(64, 256, 1024, 4096)
A  <- A_table(30, max(N7))
T7 <- t(sapply(N7, function(n) { EH <- 2 * harmonic(n) - 2; round(EH - sapply(D7, function(d) E_Rprime(n, d, A)), 3) }))
T7 <- data.frame(n = as.character(N7), T7, check.names = FALSE)
names(T7)[-1] <- paste0("d=", D7)
lim <- c(1, round(2 * (D7[-1] - 1) * sapply(D7[-1], zeta) - 1, 4))
T7 <- rbind(T7, setNames(data.frame("lim E(R')", t(lim)), names(T7)))
write.csv(T7, file.path(OUT_DIR, "SuppTable7_maximal_state.csv"), row.names = FALSE)
cat("Table 7 (E(tau'_sigma), exact)\n"); print(T7, row.names = FALSE)

## ------------------------------------------------------------------ Table 8
n  <- 256L; EH <- 2 * harmonic(n) - 2
ER_extreme <- function(d) { m <- 2:n; sum(c_n(n, m) * 2^d / m^(d - 1)) / n }  # Theorem 5

some_orientation_count <- function(X) {       # observations with an empty orthant around them
  m <- nrow(X); d <- ncol(X)
  K <- matrix(0, m, m)                        # K[i, j] = orthant of observation j around observation i
  for (q in seq_len(d)) K <- K + outer(X[, q], X[, q], "<") * 2^(q - 1)
  as.numeric(sum(vapply(seq_len(m), function(i) length(unique(K[i, -i])) < 2^d, logical(1))))
}

T8 <- do.call(rbind, lapply(c(1, 2, 3, 4, 5, 6, 8, 10, 20), function(d) {
  set.seed(SEED + 8000L + d)
  S <- 0; V <- 0
  for (m in 2:n) {
    if (m <= 2^d) { B <- m; vb <- 0 } else {
      k  <- max(40L, 6000L %/% m)
      cc <- vapply(seq_len(k), function(s) some_orientation_count(matrix(runif(m * d), m, d)), numeric(1))
      B  <- mean(cc); vb <- var(cc) / k
    }
    S <- S + c_n(n, m) * B; V <- V + c_n(n, m)^2 * vb
  }
  e <- EH - S / n; ER <- ER_extreme(d)
  data.frame(d = d, E_tau_U = round(e, 3), se = round(sqrt(V) / n, 3),
             E_tau_minus_tau_U = round(EH - ER - e, 3), E_R = round(ER, 3))
}))
write.csv(T8, file.path(OUT_DIR, "SuppTable8_maximal_state.csv"), row.names = FALSE)
cat("\nTable 8\n"); print(T8, row.names = FALSE)

## ------------------------------------------------------------------ Table 9 (tree simulation)
one_tree <- function(X) {                     # tau' and R' of the sigma-maximal state in one tree
  nn <- nrow(X); d <- ncol(X)
  H <- integer(nn); Rp <- integer(nn)
  stack <- list(list(idx = seq_len(nn), dep = 0L))
  while (length(stack)) {
    nd <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    idx <- nd$idx; m <- length(idx)
    if (m == 1L) { H[idx] <- nd$dep; next }
    P <- X[idx, , drop = FALSE]
    dom <- matrix(TRUE, m, m)                 # dom[i, j]: observation j exceeds i in every variable
    for (q in seq_len(d)) dom <- dom & outer(P[, q], P[, q], "<")
    mx <- rowSums(dom) == 0
    Rp[idx[mx]] <- Rp[idx[mx]] + 1L           # sigma-maximal in this internal node
    v <- sample.int(d, 1L); cut <- runif(1, min(P[, v]), max(P[, v]))
    left <- P[, v] < cut
    stack[[length(stack) + 1L]] <- list(idx = idx[left],  dep = nd$dep + 1L)
    stack[[length(stack) + 1L]] <- list(idx = idx[!left], dep = nd$dep + 1L)
  }
  list(tau = H - Rp, R = Rp)
}
one_sample <- function(d, s) {                # one sample of n observations, 100 trees
  set.seed(SEED + 9000L + 100L * d + s)
  X <- matrix(runif(n * d), n, d)
  tr <- lapply(seq_len(100), function(t) one_tree(X))
  Tt <- unlist(lapply(tr, `[[`, "tau")); Rr <- unlist(lapply(tr, `[[`, "R"))
  pv <- function(x) mean((x - mean(x))^2)     # variance over observations and trees
  c(mean(Tt), mean(Rr), pv(Tt), pv(Rr), 2 * mean((Tt - mean(Tt)) * (Rr - mean(Rr))))
}
run_samples <- function(d) {                  # N9 samples for one d, in parallel where possible
  if (CORES > 1L && .Platform$OS.type == "windows") {
    cl <- makeCluster(CORES); on.exit(stopCluster(cl))
    clusterExport(cl, c("one_tree", "one_sample", "n", "SEED"), envir = globalenv())
    parLapply(cl, seq_len(N9), function(s, d) one_sample(d, s), d = d)
  } else if (CORES > 1L) {
    mclapply(seq_len(N9), function(s) one_sample(d, s), mc.cores = CORES)
  } else lapply(seq_len(N9), function(s) one_sample(d, s))
}
T9 <- do.call(rbind, lapply(c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20), function(d) {
  res <- run_samples(d)
  M  <- do.call(rbind, res); mu <- colMeans(M); se <- apply(M, 2, sd) / sqrt(N9)
  cat(sprintf("d = %2d: %s  sum %.3f\n", d,
              paste(sprintf("%.3f (%.3f)", mu, se), collapse = "  "), sum(mu[3:5])))
  data.frame(d = d, E_Rprime_exact = round(E_Rprime(n, d, A), 3),
             E_tau = round(mu[1], 3), se_E_tau = round(se[1], 3),
             E_R = round(mu[2], 3), se_E_R = round(se[2], 3),
             Var_tau = round(mu[3], 3), se_Var_tau = round(se[3], 3),
             Var_R = round(mu[4], 3), se_Var_R = round(se[4], 3),
             Cov2 = round(mu[5], 3), se_Cov2 = round(se[5], 3),
             sum = round(sum(mu[3:5]), 3))
}))
write.csv(T9, file.path(OUT_DIR, "SuppTable9_maximal_state.csv"), row.names = FALSE)
cat("written to", OUT_DIR, "\n")
