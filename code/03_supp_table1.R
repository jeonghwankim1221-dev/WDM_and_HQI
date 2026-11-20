###############################################################################
# Supplementary Table 1 (simulation against exact values) from the simulation cache
#   input : OUT_DIR/sim_ch6.rds   (written by 01_exact_simul_figures.R; n = 256, 512;
#                                  d = 1..10, 20, 50, 100; 500 trees x 20 replicates)
#   output: OUT_DIR/S1_table.csv  and a short summary on the console
###############################################################################

.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))
SIM_FILE <- if (file.exists(file.path(OUT_DIR, "sim_ch6.rds"))) file.path(OUT_DIR, "sim_ch6.rds") else file.path(DATA_DIR, "sim_ch6.rds")
SIM <- readRDS(SIM_FILE)

bern_sum <- function(ps) {                     # law of a sum of independent Bernoulli(ps)
  p <- 1
  for (q in ps) p <- c(p * (1 - q), 0) + c(0, p * q)
  p
}
ext_dists <- function(n)                        # Theorem 1: H_m^ext = sum_{j<m} Bern(1/j), r = 0..m-1
  lapply(seq_len(n), function(m) if (m == 1) 1 else bern_sum(1 / seq_len(m - 1)))

kernels <- function(n, d) {                     # rows m, cols j (1..n); Theorem 4 / limit of 4.6
  Kabs <- matrix(0, n, n); KG <- matrix(0, n, n)
  for (m in 3:n) {
    j <- seq_len(m - 1)
    if (is.infinite(d)) {
      Kabs[m, 1:2] <- 2 * (1:2) / (m * (m - 1))
      if (m > 3) KG[m, 3:(m - 1)] <- 2 * (3:(m - 1)) / (m * (m - 1))
    } else {
      pm <- (2 / m)^d; am <- (2 / m)^(d - 1); den <- m * (m - 1) * (1 - pm)
      Kabs[m, 1] <- 2 * (1 - am) / den
      Kabs[m, 2:(m - 1)] <- 2 * (2 * (2 / j[-1])^(d - 1) - am) / den
      if (m > 3) KG[m, 3:(m - 1)] <- 2 * (j[3:(m - 1)] - 2 * (2 / j[3:(m - 1)])^(d - 1)) / den
    }
  }
  list(Kabs = Kabs, KG = KG)
}

joint_tau_M <- function(n, d) {                 # J[t+1, m] = P(tau = t, M_tau = m)   (Corollary 3)
  K <- kernels(n, d)
  pn <- if (is.infinite(d)) 0 else (2 / n)^d
  J <- matrix(0, n + 1, n); f <- numeric(n); J[1, n] <- pn; f[n] <- 1 - pn
  h <- 0
  while (sum(f) > 1e-16 && h < n) {
    J[h + 2, ] <- J[h + 2, ] + as.vector(f %*% K$Kabs)
    f <- as.vector(f %*% K$KG); h <- h + 1
  }
  J[seq_len(max(which(rowSums(J) > 0))), , drop = FALSE]
}

joint_tau_R <- function(J, ext) {               # P[t+1, r+1] = P(tau = t, R = r)  (Proposition 2)
  n <- ncol(J); P <- matrix(0, nrow(J), n)
  for (m in seq_len(n)) {
    w <- J[, m]; if (!any(w > 0)) next
    e <- ext[[m]]; P[, seq_along(e)] <- P[, seq_along(e)] + outer(w, e)
  }
  P
}

exact_summary <- function(n, d, ext) {
  J <- joint_tau_M(n, d); P <- joint_tau_R(J, ext)
  t <- 0:(nrow(P) - 1); r <- 0:(ncol(P) - 1)
  pt <- rowSums(P); pr <- colSums(P); m <- seq_len(n); pm <- colSums(J)
  Et <- sum(t * pt); ER <- sum(r * pr)
  Vt <- sum(t^2 * pt) - Et^2; VR <- sum(r^2 * pr) - ER^2
  Cv <- sum(outer(t, r) * P) - Et * ER
  hm <- sapply(m, function(k) sum(1 / seq_len(k - 1)))
  vm <- sapply(m, function(k) sum(1 / seq_len(k - 1)) - sum(1 / seq_len(k - 1)^2))
  hm[1] <- 0; vm[1] <- 0
  c(n = n, d = d, E_tau = Et, E_R = ER, E_H = Et + ER, Var_tau = Vt, Var_R = VR, Cov = Cv,
    Var_H = Vt + VR + 2 * Cv, E_M = sum(m * pm), E_v = sum(vm * pm),
    Var_h = sum(hm^2 * pm) - sum(hm * pm)^2)
}


keys <- c("E_tau", "E_R", "Var_tau", "Var_R", "Cov")
rows <- list()
for (nm in names(SIM)) {
  nd <- as.integer(strsplit(nm, " ")[[1]]); n <- nd[1]; d <- nd[2]
  ex <- exact_summary(n, d, ext_dists(n))
  m <- SIM[[nm]]$mean; se <- SIM[[nm]]$se
  z <- sapply(keys, function(k) if (se[[k]] > 0) (m[[k]] - ex[[k]]) / se[[k]] else NA)
  rows[[nm]] <- data.frame(n = n, d = d, t(round(z, 2)),
                           dE_H = round(m[["E_H"]] - ex[["E_H"]], 3), dVar_H = round(m[["Var_H"]] - ex[["Var_H"]], 3),
                           dE_tau = m[["E_tau"]] - ex[["E_tau"]], dE_R = m[["E_R"]] - ex[["E_R"]],
                           dVar_R = m[["Var_R"]] - ex[["Var_R"]],
                           dE_v = m[["E_v"]] - ex[["E_v"]], se_E_v = se[["E_v"]],
                           dVar_h = m[["Var_h"]] - ex[["Var_h"]], se_Var_h = se[["Var_h"]])
}
TAB <- do.call(rbind, rows); TAB <- TAB[order(TAB$n, TAB$d), ]
write.csv(TAB, file.path(OUT_DIR, "SuppTable1_simulation_vs_exact.csv"), row.names = FALSE)
print(TAB[, c("n", "d", keys, "dE_H", "dVar_H")], row.names = FALSE)
cat(sprintf("\nsettings: %d, comparisons (non-missing z): %d\n", nrow(TAB), sum(!is.na(as.matrix(TAB[, keys])))))
cat(sprintf("max |z|                      : %.2f\n", max(abs(as.matrix(TAB[, keys])), na.rm = TRUE)))
cat(sprintf("max |diff| E(tau), E(R)      : %.3f\n", max(abs(c(TAB$dE_tau, TAB$dE_R)))))
cat(sprintf("max |diff| Var(R)            : %.3f\n", max(abs(TAB$dVar_R))))
cat(sprintf("max |diff| E(H), Var(H)      : %.3f, %.3f\n", max(abs(TAB$dE_H)), max(abs(TAB$dVar_H))))
cat(sprintf("d = 1: E{v} diff %.3f (se %.3f); Var{h} diff %.3f (se %.3f)  [n = 256]\n",
            TAB$dE_v[TAB$n == 256 & TAB$d == 1], TAB$se_E_v[TAB$n == 256 & TAB$d == 1],
            TAB$dVar_h[TAB$n == 256 & TAB$d == 1], TAB$se_Var_h[TAB$n == 256 & TAB$d == 1]))
