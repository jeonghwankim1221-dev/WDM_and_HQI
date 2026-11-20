###############################################################################
# Core functions for the real-data analysis (Results: real data; Table 2 and Extended Data Table 2)
# Extracted from IF_realdata_components_v3.R and parameterised so that
# 04_realdata.R can call it once per (dataset, marginal) pair.
###############################################################################
Hn  <- function(n) sum(1 / seq_len(n))
Hn2 <- function(n) sum(1 / seq_len(n)^2)

#  tree: H, tau, M_tau for every observation (splitting rule of Methods) ------
one_tree <- function(X) {
  n <- nrow(X); d <- ncol(X)
  H <- integer(n); TAU <- rep(NA_integer_, n); MT <- rep(NA_integer_, n)
  stack <- vector("list", 2L * n); top <- 1L; stack[[1L]] <- list(idx = seq_len(n), h = 0L)
  while (top > 0L) {
    nd <- stack[[top]]; top <- top - 1L; idx <- nd$idx; h <- nd$h; m <- length(idx)
    if (m == 1L) { H[idx] <- h; if (is.na(TAU[idx])) { TAU[idx] <- h; MT[idx] <- 1L }; next }
    S <- X[idx, , drop = FALSE]
    v1 <- S[, 1L]; cand <- idx[v1 == min(v1) | v1 == max(v1)]; cand <- cand[is.na(TAU[cand])]
    for (cc in cand) { row <- rep(X[cc, ], each = m)
      if (all(colSums(S < row) == 0L | colSums(S > row) == 0L)) { TAU[cc] <- h; MT[cc] <- m } }
    q <- sample.int(d, 1L); v <- S[, q]; rg <- range(v); cut <- runif(1L, rg[1L], rg[2L])
    L <- idx[v < cut]; R <- idx[v >= cut]
    if (length(L) == 0L || length(R) == 0L) { top <- top + 1L; stack[[top]] <- nd; next }
    top <- top + 1L; stack[[top]] <- list(idx = L, h = h + 1L)
    top <- top + 1L; stack[[top]] <- list(idx = R, h = h + 1L)
  }
  list(H = H, tau = TAU, M = MT)
}

# theory for (n, d): recursion of Corollary 1 -------------------------------
theory <- function(n, d) {
  pm <- function(m) (2 / m)^d; am <- function(m) (2 / m)^(d - 1)
  P <- matrix(0, n + 1, n + 1)
  f <- numeric(n + 1); f[n + 1] <- 1 - pm(n); P[n + 1, 1] <- pm(n); h <- 0L
  while (sum(f) > 1e-15 && h < n) {
    g <- numeric(n + 1)
    for (m in which(f > 0) - 1L) { if (m < 3) next
      j <- 1:(m - 1); den <- m * (m - 1) * (1 - pm(m))
      K <- 2 * (j - am(m)) / den; Kab <- 2 * (2 * (2 / j)^(d - 1) - am(m)) / den; Kab[1] <- K[1]
      KG <- K - Kab; KG[j <= 2] <- 0
      P[j + 1, h + 2] <- P[j + 1, h + 2] + f[m + 1] * Kab; g[j + 1] <- g[j + 1] + f[m + 1] * KG }
    f <- g; h <- h + 1L }
  m <- 0:n; t <- 0:n
  hM <- sapply(m, function(k) if (k <= 1) 0 else Hn(k - 1)); vM <- sapply(m, function(k) if (k <= 1) 0 else Hn(k - 1) - Hn2(k - 1))
  pM <- rowSums(P); pt <- colSums(P)
  Et <- sum(pt * t); Vt <- sum(pt * t^2) - Et^2; ER <- sum(pM * hM); Ev <- sum(pM * vM); Vh <- sum(pM * hM^2) - ER^2
  cov <- sum(outer(hM, t) * P) - Et * ER
  c(E_tau = Et, E_R = ER, Var_tau = Vt, Ev = Ev, Vh = Vh, Cov = cov)
}
limit_vals <- function(n) { VarH <- 2 * Hn(n) + 2 - 4 * Hn2(n)
  c(E_tau = 2 * Hn(n) - 2 - 2 / 3, Var_tau = VarH - 2 / 9) }
Dcomp <- function(v, lim) c(D_tau = unname(v["Var_tau"] - lim["Var_tau"]), D_v = unname(v["Ev"]),
                            D_h = unname(v["Vh"] - 2 / 9), D_c = unname(2 * v["Cov"]), delta = unname(v["E_tau"] - lim["E_tau"]))

# empirical components from one subsample of X_all (psi drawn without replacement)
one_sub <- function(seed, X_all, N, psi, n_trees) {
  set.seed(seed); X <- X_all[sample.int(N, psi), , drop = FALSE]
  tau <- H <- M <- integer(0)
  for (t in seq_len(n_trees)) { r <- one_tree(X); tau <- c(tau, r$tau); H <- c(H, r$H); M <- c(M, r$M) }
  R <- H - tau
  grp <- split(R, M); w <- sapply(grp, length) / length(R)
  Ev <- sum(w * sapply(grp, function(x) mean((x - mean(x))^2)))
  Vh <- sum(w * (sapply(grp, mean) - mean(R))^2)
  c(E_tau = mean(tau), E_R = mean(R), Var_tau = mean((tau - mean(tau))^2), Ev = Ev, Vh = Vh,
    Cov = mean((tau - mean(tau)) * (R - mean(R))))
}

# Read a data set from /data (no download: the capsule ships the five CSV files) -
load_dataset_local <- function(name, data_dir) {
  csv <- file.path(data_dir, paste0(name, ".csv"))
  if (!file.exists(csv)) stop("missing data file: ", csv,
    "\nPlace quakes.csv, magic.csv, wine_white.csv, letter.csv, supercond.csv in /data (explanatory variables only).")
  X <- as.matrix(read.csv(csv)); X <- X[, setdiff(colnames(X), "y"), drop = FALSE]
  storage.mode(X) <- "double"; X
}

# Full pipeline for one (dataset, marginal): jitter -> (rank transform) -> psi subsamples -> trees
run_one_dataset <- function(name, marginal, data_dir, psi, n_sub, n_trees, jitter_sd, seed, cores) {
  X_all <- load_dataset_local(name, data_dir)
  N <- nrow(X_all); d <- ncol(X_all)
  set.seed(seed)
  sds <- apply(X_all, 2, sd); sds[!is.finite(sds) | sds == 0] <- 1
  X_all <- X_all + matrix(rnorm(N * d, 0, jitter_sd), N, d) * rep(sds, each = N)         # tie breaking
  if (marginal == "rank") X_all <- apply(X_all, 2, function(x) rank(x) / (N + 1))        # uniform marginals

  seeds <- seed + seq_len(n_sub)
  res <- if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("one_tree", "one_sub", "X_all", "N", "psi", "n_trees"), envir = environment())
    r <- parallel::parLapply(cl, seeds, function(s) one_sub(s, X_all, N, psi, n_trees)); parallel::stopCluster(cl); r
  } else parallel::mclapply(seeds, function(s) one_sub(s, X_all, N, psi, n_trees), mc.cores = cores)
  EMP <- do.call(rbind, res); emp <- colMeans(EMP); se <- apply(EMP, 2, sd) / sqrt(nrow(EMP))

  TH <- theory(psi, d); LIM <- limit_vals(psi)
  D_emp <- Dcomp(emp, LIM); D_th <- Dcomp(TH, LIM)
  D_se  <- c(D_tau = se[["Var_tau"]], D_v = se[["Ev"]], D_h = se[["Vh"]], D_c = 2 * se[["Cov"]], delta = se[["E_tau"]])

  mom <- data.frame(name = names(TH), theory = unname(TH), empirical = unname(emp), se = unname(se))
  Dtab <- data.frame(name = names(D_th), theory = unname(D_th), empirical = unname(D_emp), se = unname(D_se))
  list(dataset = name, marginal = marginal, N = N, d = d,
       moment = mom, D = Dtab)
}
