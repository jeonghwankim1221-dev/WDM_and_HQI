###############################################################################
# Post-absorption law under rank dependence — v3
#
# v3 fixes an aggregation bug in v2: aggregate(cbind(...) ~ .) includes the
# left-hand variables in the dot expansion, which silently replaced sums by
# group values. All sum/count aggregation now uses rowsum() on explicit keys.
#
# Changes from v1
#   - n = 256, 300 trees, 40 replications: about 25x more entry records at
#     m = 8 than v1 (node counts scale with n, c_n(m) = 2n/(m(m+1))).
#   - Kernels and TV are computed from counts POOLED across replications.
#     Per-replication TV in v1 was dominated by multinomial sampling noise.
#   - Closure test by parent size: for each special record we store the size
#     of the parent node. If the size-only recursion closes, the residual
#     from a size-k special node must not depend on the parent's size:
#         r_{k|p} = E[Rrem | size k, parent size p] = r_k  for all p.
#     This is the direct test of H5-2 (one-step Markov in size); by the
#     Markov property it covers longer histories as well.
#   - Theory reference for E[H] computed per n (v1 used the first n only).
#
# Copulas (exchangeable d-variate), Kendall tau_K = 0.5:
#   indep, gaussian (lambda = sin(pi/4)), t3 (same lambda, nu = 3),
#   clayton (theta = 2), sclayton (survival), gumbel (theta = 2), sgumbel.
#
# Records at every internal node A, for every point i extremal in A:
#   m = |A|, entry = (tau_i == depth A), J = child size holding i,
#   Rrem = H_i - depth(A), tau_i, psize = |parent(A)| (NA at root).
#
# R 4.3.3, base R + parallel, Linux/Windows. RDS-first output next to script.
###############################################################################

rm(list = ls())
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))
suppressWarnings(suppressMessages(library(parallel)))

CFG <- list(
  copulas  = c("indep", "gaussian", "t3", "clayton", "sclayton", "gumbel", "sgumbel"),
  tau_K    = 0.5,
  dims     = c(2L, 3L),
  n        = 256L,
  n_trees  = 300L,
  n_rep    = 40L,
  m_focus  = c(4L, 6L, 8L, 12L, 16L),
  psize_bins = c(0, 4, 8, 16, 32, 64, 128, 256, Inf),
  n_cores  = max(1L, parallel::detectCores() - 1L),
  seed     = 20260917L
)
LAMBDA  <- sin(pi * CFG$tau_K / 2)
THETA_C <- 2 * CFG$tau_K / (1 - CFG$tau_K)
THETA_G <- 1 / (1 - CFG$tau_K)

RNGkind("L'Ecuyer-CMRG"); set.seed(CFG$seed)

par_lapply <- function(X, FUN, cores, export, seed) {
  if (cores <= 1L) return(lapply(X, FUN))
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores); on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, seed)
    parallel::clusterExport(cl, export, envir = .GlobalEnv)
    parallel::parLapplyLB(cl, X, FUN)
  } else parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
}
pick_one <- function(v) v[sample.int(length(v), 1L)]
se <- function(x) { x <- x[is.finite(x)]; if (length(x) <= 1L) NA_real_ else sd(x) / sqrt(length(x)) }
harm <- function(N) c(0, cumsum(1 / seq_len(N)))

## ---------------------------------------------------------------------------
## copula samplers (uniform marginals, exchangeable)
## ---------------------------------------------------------------------------
rpos_stable <- function(n, alpha) {
  W <- runif(n, 0, pi); E <- rexp(n)
  (sin(alpha * W) / sin(W)^(1 / alpha)) * (sin((1 - alpha) * W) / E)^((1 - alpha) / alpha)
}
rcop <- function(fam, n, d) {
  switch(fam,
    indep    = matrix(runif(n * d), n, d),
    gaussian = { z0 <- rnorm(n)
                 pnorm(sqrt(LAMBDA) * z0 + sqrt(1 - LAMBDA) * matrix(rnorm(n * d), n, d)) },
    t3       = { nu <- 3; z0 <- rnorm(n)
                 Z <- sqrt(LAMBDA) * z0 + sqrt(1 - LAMBDA) * matrix(rnorm(n * d), n, d)
                 pt(Z * sqrt(nu / rchisq(n, nu)), df = nu) },
    clayton  = { V <- rgamma(n, shape = 1 / THETA_C, rate = 1)
                 (1 + matrix(rexp(n * d), n, d) / V)^(-1 / THETA_C) },
    sclayton = 1 - rcop("clayton", n, d),
    gumbel   = { V <- rpos_stable(n, 1 / THETA_G)
                 exp(-(matrix(rexp(n * d), n, d) / V)^(1 / THETA_G)) },
    sgumbel  = 1 - rcop("gumbel", n, d),
    stop("unknown copula ", fam))
}

## ---------------------------------------------------------------------------
## tree with special-point records (incl. parent size)
## ---------------------------------------------------------------------------
build_tree <- function(X) {
  n <- nrow(X); d <- ncol(X)
  ord <- seq_len(n); H <- integer(n); tau <- rep(NA_integer_, n)
  sb <- integer(2L * n); sE <- integer(2L * n); sD <- integer(2L * n); sP <- integer(2L * n)
  top <- 1L; sb[1] <- 1L; sE[1] <- n; sD[1] <- 0L; sP[1] <- NA_integer_
  cap <- 4L * n
  Rm <- integer(cap); Re <- logical(cap); RJ <- integer(cap)
  Ri <- integer(cap); Rd <- integer(cap); Rp <- integer(cap); nr <- 0L
  viol <- 0L
  while (top > 0L) {
    b <- sb[top]; e <- sE[top]; dep <- sD[top]; psz <- sP[top]; top <- top - 1L
    m <- e - b + 1L; idx <- ord[b:e]
    if (m == 1L) { if (is.na(tau[idx])) tau[idx] <- dep; H[idx] <- dep; next }
    mins <- numeric(d); maxs <- numeric(d); ext <- rep(TRUE, m)
    for (q in seq_len(d)) {
      xq <- X[idx, q]; mn <- min(xq); mx <- max(xq); mins[q] <- mn; maxs[q] <- mx
      if (m > 2L) ext <- ext & (xq == mn | xq == mx)
    }
    prev <- !is.na(tau[idx]) & (tau[idx] < dep)
    if (any(prev)) viol <- viol + sum(prev & !ext)
    newly <- ext & is.na(tau[idx]); tau[idx[newly]] <- dep
    valid <- which(maxs > mins)
    if (!length(valid)) { H[idx] <- dep; next }
    q <- pick_one(valid); s <- runif(1L, mins[q], maxs[q]); gl <- X[idx, q] <= s
    if (all(gl) || !any(gl)) {
      ux <- sort(unique(X[idx, q])); k <- pick_one(seq_len(length(ux) - 1L))
      gl <- X[idx, q] <= (ux[k] + ux[k + 1L]) / 2
    }
    nL <- sum(gl)
    for (j in which(ext)) {
      nr <- nr + 1L
      Rm[nr] <- m; Re[nr] <- newly[j]; RJ[nr] <- if (gl[j]) nL else m - nL
      Ri[nr] <- idx[j]; Rd[nr] <- dep; Rp[nr] <- psz
    }
    ord[b:e] <- c(idx[gl], idx[!gl])
    top <- top + 1L; sb[top] <- b;      sE[top] <- b + nL - 1L; sD[top] <- dep + 1L; sP[top] <- m
    top <- top + 1L; sb[top] <- b + nL; sE[top] <- e;           sD[top] <- dep + 1L; sP[top] <- m
  }
  k <- seq_len(nr)
  list(H = H, tau = tau, viol = viol,
       rec = data.frame(m = Rm[k], entry = Re[k], J = RJ[k], Rrem = H[Ri[k]] - Rd[k],
                        tau_i = tau[Ri[k]], psize = Rp[k]))
}

## ---------------------------------------------------------------------------
## one cell: aggregate to sums/counts only
## ---------------------------------------------------------------------------
run_cell <- function(fam, d, n, r) {
  X <- rcop(fam, n, d)
  recs <- vector("list", CFG$n_trees); viol <- 0L
  Hs <- numeric(CFG$n_trees * n); ts <- numeric(CFG$n_trees * n)
  for (t in seq_len(CFG$n_trees)) {
    tr <- build_tree(X); recs[[t]] <- tr$rec; viol <- viol + tr$viol
    sl <- ((t - 1L) * n + 1L):(t * n); Hs[sl] <- tr$H; ts[sl] <- tr$tau
  }
  rec <- do.call(rbind, recs)
  rec$pbin <- cut(rec$psize, CFG$psize_bins, labels = FALSE)
  rec$pbin[is.na(rec$psize)] <- 0L

  # explicit-key aggregation (rowsum), immune to formula dot-expansion
  agg3 <- function(df, keys) {
    key <- do.call(paste, c(df[keys], sep = "|"))
    M <- rowsum(cbind(s1 = df$Rrem, s2 = df$Rrem^2, cnt = rep(1, nrow(df))), key)
    parts <- do.call(rbind, strsplit(rownames(M), "|", fixed = TRUE))
    out <- as.data.frame(M)
    for (i in seq_along(keys)) out[[keys[i]]] <- type.convert(parts[, i], as.is = TRUE)
    rownames(out) <- NULL
    out[c(keys, "s1", "s2", "cnt")]
  }
  aggc <- function(df, keys) {
    key <- do.call(paste, c(df[keys], sep = "|"))
    M <- rowsum(cbind(cnt = rep(1, nrow(df))), key)
    parts <- do.call(rbind, strsplit(rownames(M), "|", fixed = TRUE))
    out <- as.data.frame(M)
    for (i in seq_along(keys)) out[[keys[i]]] <- type.convert(parts[, i], as.is = TRUE)
    rownames(out) <- NULL
    out[c(keys, "cnt")]
  }
  by_m     <- agg3(rec, c("m", "entry"))
  by_kern  <- aggc(rec, c("m", "J", "entry"))
  by_hist  <- agg3(rec[rec$entry, ], c("m", "tau_i"))
  by_par   <- agg3(rec, c("m", "pbin"))
  tag <- function(df) { df$copula <- fam; df$d <- d; df$n <- n; df$rep <- r; df }
  list(by_m = tag(by_m), by_kern = tag(by_kern), by_hist = tag(by_hist), by_par = tag(by_par),
       glob = tag(data.frame(mean_H = mean(Hs), mean_tau = mean(ts), mean_R = mean(Hs - ts),
                             viol = viol)))
}

## ---------------------------------------------------------------------------
## run
## ---------------------------------------------------------------------------
conc <- do.call(rbind, lapply(CFG$copulas, function(f) {
  U <- rcop(f, 30000L, 2L)
  data.frame(copula = f, kendall = cor(U[, 1], U[, 2], method = "kendall"),
             spearman = cor(U[, 1], U[, 2], method = "spearman"),
             lower_tail_0.02 = mean(U[, 1] < 0.02 & U[, 2] < 0.02) / 0.02,
             upper_tail_0.02 = mean(U[, 1] > 0.98 & U[, 2] > 0.98) / 0.02)
}))

jobs <- list()
for (f in CFG$copulas) for (d in CFG$dims) for (r in seq_len(CFG$n_rep))
  jobs[[length(jobs) + 1L]] <- list(f = f, d = d, n = CFG$n, r = r)
message(sprintf("%d cells, %d cores", length(jobs), CFG$n_cores))
EXPORT <- c("CFG", "LAMBDA", "THETA_C", "THETA_G", "pick_one", "se", "harm",
            "rpos_stable", "rcop", "build_tree", "run_cell", "jobs")
t0 <- Sys.time()
res <- par_lapply(seq_along(jobs), function(i) { jb <- jobs[[i]]; run_cell(jb$f, jb$d, jb$n, jb$r) },
                  cores = CFG$n_cores, export = EXPORT, seed = CFG$seed)
bad <- vapply(res, inherits, logical(1), "try-error")
if (any(bad)) stop("cell failed: ", res[[which(bad)[1]]])
bind <- function(nm) do.call(rbind, lapply(res, `[[`, nm))
BYM <- bind("by_m"); KERN <- bind("by_kern"); HIST <- bind("by_hist"); PAR <- bind("by_par"); GLOB <- bind("glob")

## ---------------------------------------------------------------------------
## summaries
## ---------------------------------------------------------------------------
Hn <- harm(CFG$n + 2L)

# residual means: pooled over reps, plus rep-level SE
pool_mean <- function(df, keys) {
  key <- do.call(paste, c(df[keys], sep = "|"))
  M <- rowsum(as.matrix(df[c("s1", "s2", "cnt")]), key)
  parts <- do.call(rbind, strsplit(rownames(M), "|", fixed = TRUE))
  p <- as.data.frame(M)
  for (i in seq_along(keys)) p[[keys[i]]] <- type.convert(parts[, i], as.is = TRUE)
  rownames(p) <- NULL
  p <- p[c(keys, "s1", "s2", "cnt")]
  p$mean <- p$s1 / p$cnt
  p$sd_within <- sqrt(pmax(p$s2 / p$cnt - p$mean^2, 0))
  p
}
rep_se <- function(df, keys) {
  df$mu <- df$s1 / df$cnt
  do.call(rbind, lapply(split(df, df[keys], drop = TRUE), function(g)
    cbind(g[1, keys, drop = FALSE], n_rep = nrow(g), rep_se = se(g$mu), mean_count = mean(g$cnt))))
}
R_POOL <- pool_mean(BYM, c("copula", "d", "n", "m", "entry"))
R_POOL <- merge(R_POOL, rep_se(BYM, c("copula", "d", "n", "m", "entry")), by = c("copula", "d", "n", "m", "entry"))
R_POOL$H_ref <- Hn[R_POOL$m]
R_POOL$Delta <- R_POOL$mean - R_POOL$H_ref
R_ENT <- R_POOL[R_POOL$entry, ]
# all-node residual (entry + descendant) by m
R_ALL <- pool_mean(BYM, c("copula", "d", "n", "m"))
R_ALL$H_ref <- Hn[R_ALL$m]; R_ALL$Delta <- R_ALL$mean - R_ALL$H_ref

# kernels pooled over reps
K_POOL <- { key <- do.call(paste, c(KERN[c("copula", "d", "n", "m", "J", "entry")], sep = "|"))
            M <- rowsum(cbind(cnt = KERN$cnt), key)
            parts <- do.call(rbind, strsplit(rownames(M), "|", fixed = TRUE))
            k <- data.frame(copula = parts[, 1], d = as.integer(parts[, 2]), n = as.integer(parts[, 3]),
                            m = as.integer(parts[, 4]), J = as.integer(parts[, 5]),
                            entry = as.logical(parts[, 6]), cnt = as.numeric(M[, 1]))
            k }
kernel_shape <- function(kt) {
  do.call(rbind, lapply(split(kt, kt[c("copula", "d", "n", "m", "entry")], drop = TRUE), function(g) {
    m <- g$m[1]; if (m < 3L) return(NULL)
    full <- numeric(m - 1L); full[g$J] <- g$cnt; tot <- sum(full); p <- full / tot
    j <- seq_len(m - 1L)
    data.frame(copula = g$copula[1], d = g$d[1], n = g$n[1], m = m, entry = g$entry[1],
               total = tot, TV = 0.5 * sum(abs(p - 1 / (m - 1))),
               TV_noise_floor = 0.5 * (m - 1) * sqrt((1 / (m - 1)) * (1 - 1 / (m - 1)) / tot) * sqrt(2 / pi),
               EJ = sum(j * p), EJ_unif = m / 2,
               VarJ = sum(j^2 * p) - sum(j * p)^2, VarJ_unif = ((m - 1)^2 - 1) / 12,
               p_J1 = p[1], p_J1_unif = 1 / (m - 1))
  }))
}
SHAPE <- kernel_shape(K_POOL)

# closure Gamma with pooled all-node kernel and pooled all-node residuals
rmap_of <- function(cop, d, n) { s <- R_ALL[R_ALL$copula == cop & R_ALL$d == d & R_ALL$n == n, ]; setNames(s$mean, s$m) }
GAMMA <- do.call(rbind, lapply(split(K_POOL, K_POOL[c("copula", "d", "n", "m")], drop = TRUE), function(g) {
  m <- g$m[1]; if (m < 3L) return(NULL)
  ga <- aggregate(cnt ~ J, data = g, FUN = sum); p <- ga$cnt / sum(ga$cnt)
  rm <- rmap_of(g$copula[1], g$d[1], g$n[1])
  rk <- rm[as.character(ga$J)]; rk[is.na(rk)] <- 0; rk[ga$J == 1L] <- 0
  obs <- rm[as.character(m)]
  data.frame(copula = g$copula[1], d = g$d[1], n = g$n[1], m = m,
             r_obs = obs, r_pred = 1 + sum(p * rk), Gamma = obs - (1 + sum(p * rk)))
}))
# rep-level SE for Gamma from per-rep residual means (approximate: kernel fixed at pooled)
GAMMA_SE <- do.call(rbind, lapply(split(BYM, BYM[c("copula", "d", "n", "rep")], drop = TRUE), function(g) {
  A <- rowsum(cbind(g$s1, g$cnt), g$m); rm <- setNames(A[, 1] / A[, 2], rownames(A))
  kp <- K_POOL[K_POOL$copula == g$copula[1] & K_POOL$d == g$d[1] & K_POOL$n == g$n[1], ]
  do.call(rbind, lapply(split(kp, kp$m), function(k) {
    m <- k$m[1]; if (m < 3L || is.na(rm[as.character(m)])) return(NULL)
    ga <- aggregate(cnt ~ J, data = k, FUN = sum); p <- ga$cnt / sum(ga$cnt)
    rk <- rm[as.character(ga$J)]; rk[is.na(rk)] <- 0; rk[ga$J == 1L] <- 0
    data.frame(copula = g$copula[1], d = g$d[1], n = g$n[1], rep = g$rep[1], m = m,
               Gamma_rep = rm[as.character(m)] - (1 + sum(p * rk)))
  }))
}))
GS <- do.call(rbind, lapply(split(GAMMA_SE, GAMMA_SE[c("copula", "d", "n", "m")], drop = TRUE), function(g)
  data.frame(copula = g$copula[1], d = g$d[1], n = g$n[1], m = g$m[1], Gamma_se = se(g$Gamma_rep))))
GAMMA <- merge(GAMMA, GS, by = c("copula", "d", "n", "m"))
GAMMA$z <- GAMMA$Gamma / GAMMA$Gamma_se

# parent-size closure test: r_{k|p} vs r_k
PARP <- pool_mean(PAR, c("copula", "d", "n", "m", "pbin"))
PARP <- merge(PARP, R_ALL[, c("copula", "d", "n", "m", "mean")], by = c("copula", "d", "n", "m"), suffixes = c("", "_all"))
PARP$diff <- PARP$mean - PARP$mean_all
PARP$se_within <- PARP$sd_within / sqrt(PARP$cnt)
PARP$z_within <- PARP$diff / PARP$se_within
PBIN_LABEL <- c("root", "2-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256", ">256")
PARP$parent_bin <- PBIN_LABEL[PARP$pbin + 1L]
PARSUM <- do.call(rbind, lapply(split(PARP, PARP[c("copula", "d", "n", "m")], drop = TRUE), function(g) {
  g <- g[g$cnt >= 200, ]
  if (nrow(g) < 2L) return(NULL)
  data.frame(copula = g$copula[1], d = g$d[1], n = g$n[1], m = g$m[1], n_bins = nrow(g),
             chi2 = sum((g$diff / g$se_within)^2), df = nrow(g) - 1L,
             max_abs_diff = max(abs(g$diff)), spread = max(g$mean) - min(g$mean))
}))
PARSUM$chi2_per_df <- PARSUM$chi2 / PARSUM$df

# history
HISTP <- pool_mean(HIST, c("copula", "d", "n", "m", "tau_i"))
HISTP$Delta <- HISTP$mean - Hn[HISTP$m]; HISTP$se_within <- HISTP$sd_within / sqrt(HISTP$cnt)

# global
GLOBS <- do.call(rbind, lapply(split(GLOB, GLOB[c("copula", "d", "n")], drop = TRUE), function(g)
  data.frame(copula = g$copula[1], d = g$d[1], n = g$n[1], mean_H = mean(g$mean_H), mean_H_se = se(g$mean_H),
             mean_tau = mean(g$mean_tau), mean_R = mean(g$mean_R), mean_R_se = se(g$mean_R))))
GLOBS$theory_E_H <- vapply(GLOBS$n, function(nn) 2 * harm(nn)[nn + 1L] - 2, numeric(1))
GLOBS$dev_E_H <- GLOBS$mean_H - GLOBS$theory_E_H

## ---------------------------------------------------------------------------
## output
## ---------------------------------------------------------------------------
RESULT <- list(config = CFG, params = list(lambda = LAMBDA, theta_clayton = THETA_C, theta_gumbel = THETA_G),
               concordance = conc, residual_entry = R_ENT, residual_all = R_ALL, kernel_shape = SHAPE,
               gamma = GAMMA, parent_test = PARP, parent_summary = PARSUM, history = HISTP, global = GLOBS,
               raw = list(by_m = BYM, kern = KERN, hist = HIST, par = PAR, glob = GLOB),
               viol_total = sum(GLOB$viol), elapsed_min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
saveRDS(RESULT, file.path(OUT_DIR, "copula_result.rds"))
write.csv(R_ENT, file.path(OUT_DIR, "SuppTable3_copula_entry.csv"), row.names = FALSE)
write.csv(SHAPE, file.path(OUT_DIR, "copula_kernel_shape.csv"), row.names = FALSE)
write.csv(GAMMA, file.path(OUT_DIR, "copula_gamma.csv"), row.names = FALSE)
write.csv(PARSUM, file.path(OUT_DIR, "copula_parent_summary.csv"), row.names = FALSE)
write.csv(PARP, file.path(OUT_DIR, "copula_parent_detail.csv"), row.names = FALSE)
write.csv(HISTP, file.path(OUT_DIR, "copula_history.csv"), row.names = FALSE)
write.csv(GLOBS, file.path(OUT_DIR, "copula_global.csv"), row.names = FALSE)
write.csv(conc, file.path(OUT_DIR, "copula_concordance.csv"), row.names = FALSE)

cat("\n=======================================================================\n")
cat(sprintf("elapsed %.1f min, absorption violations = %d\n", RESULT$elapsed_min, RESULT$viol_total))
cat("\nconcordance\n"); print(conc, digits = 4, row.names = FALSE)
cat("\nglobal E[H] vs 2H_n - 2\n")
print(GLOBS[order(GLOBS$d, GLOBS$copula), c("copula", "d", "mean_H", "mean_H_se", "dev_E_H", "mean_R")], digits = 5, row.names = FALSE)
show <- function(df, cols, title) {
  cat("\n", title, "\n", sep = ""); s <- df[df$m %in% CFG$m_focus, ]
  print(s[order(s$d, s$m, s$copula), cols], digits = 4, row.names = FALSE)
}
show(R_ENT, c("copula", "d", "m", "mean", "Delta", "rep_se", "cnt"), "Delta(m) = E[R | entry at m] - H_{m-1}  (pooled; rep_se across reps)")
show(GAMMA, c("copula", "d", "m", "r_obs", "r_pred", "Gamma", "Gamma_se", "z"), "closure Gamma(m), pooled kernel")
show(SHAPE[SHAPE$entry, ], c("copula", "d", "m", "total", "TV", "TV_noise_floor", "EJ", "p_J1", "p_J1_unif"),
     "entry kernel shape (pooled): TV vs its sampling floor, E[J|m] (unif m/2), P(J=1) (unif 1/(m-1))")
show(PARSUM, c("copula", "d", "m", "n_bins", "chi2_per_df", "max_abs_diff", "spread"),
     "parent-size closure: r_{k|p} vs r_k  (chi2/df ~ 1 if size-only closure holds)")
cat("\nhistory at m = 8, d = 3: Delta by entry depth\n")
hh <- HISTP[HISTP$m == 8L & HISTP$d == 3L & HISTP$cnt >= 100, ]
print(hh[order(hh$copula, hh$tau_i), c("copula", "tau_i", "Delta", "se_within", "cnt")], digits = 4, row.names = FALSE)
cat("\nwritten to:", OUT_DIR, "\n=======================================================================\n")
