# 09_supp_table6_zeta_limits.R
# Recursion of Corollary 1 at large n (default 4096) versus the n -> Inf limits of
# Supplementary Note 11 (Supplementary Table 6). Self-contained: no other script needed.
#
# For each d, the joint law of (tau, M_tau) is computed by iterating the reduced kernels
# of Theorem 4; moments of R and Cov(tau, R) follow from h(m) = H_{m-1}, v(m) = H_{m-1} - H^{(2)}_{m-1}
# (Proposition 2), so the extreme-state depth distributions are not needed.
#
# Memory: two dense n x n kernels (n = 4096 -> about 270 MB). Runtime a few minutes.

n     <- as.integer(Sys.getenv("N_CHECK", "4096"))
d_set <- 1:10
OUT_DIR <- Sys.getenv("OUT_DIR", ".")   # set by code/run (/results); current folder otherwise

harm  <- function(k, p = 1) cumsum(c(0, (1:max(k, 1))^(-p)))[k + 1]   # H_k^{(p)}, H_0 = 0
Hn    <- harm(n); Hn2 <- harm(n, 2)
E_H   <- 2 * Hn - 2
Var_H <- 2 * Hn + 2 - 4 * Hn2
m     <- 1:n
h_m   <- harm(m - 1)                      # h(m) = H_{m-1}
v_m   <- harm(m - 1) - harm(m - 1, 2)     # v(m)

## Supplementary Table 6 (limits as n -> Inf), columns: delta, D_tau, D_v, D_h, D_c, Var_R
tab6 <- rbind(
  c(-1.333333, 2.222222, 0.710132, 1.067646, -4.000000, 2.000000),
  c(-0.492806, 0.331070, 0.196191, 0.310648, -0.837909, 0.729061),
  c(-0.247299, 0.091788, 0.083786, 0.137459, -0.313033, 0.443467),
  c(-0.139746, 0.032671, 0.042926, 0.071882, -0.147479, 0.337030),
  c(-0.083884, 0.013234, 0.024194, 0.040941, -0.078368, 0.287357),
  c(-0.052143, 0.005780, 0.014421, 0.024508, -0.044709, 0.261151),
  c(-0.033128, 0.002646, 0.008902, 0.015141, -0.026689, 0.246265),
  c(-0.021350, 0.001249, 0.005623, 0.009553, -0.016425, 0.237398),
  c(-0.013894, 0.000602, 0.003608, 0.006116, -0.010326, 0.231946),
  c(-0.009102, 0.000294, 0.002340, 0.003957, -0.006591, 0.228519))
colnames(tab6) <- c("delta", "D_tau", "D_v", "D_h", "D_c", "Var_R")

kernels <- function(n, d) {
  Kabs <- matrix(0, n, n); KG <- matrix(0, n, n)        # row = parent size m, col = child size j
  for (mm in 3:n) {
    pm <- (2 / mm)^d; am <- (2 / mm)^(d - 1); den <- mm * (mm - 1) * (1 - pm)
    j <- 1:(mm - 1)
    kabs <- numeric(mm - 1)
    kabs[1] <- 2 * (1 - am) / den
    if (mm > 2) kabs[2:(mm - 1)] <- 2 * (2 * (2 / j[-1])^(d - 1) - am) / den
    kg <- numeric(mm - 1)
    if (mm > 3) kg[3:(mm - 1)] <- 2 * (j[3:(mm - 1)] - 2 * (2 / j[3:(mm - 1)])^(d - 1)) / den
    Kabs[mm, j] <- kabs; KG[mm, j] <- kg
  }
  list(Kabs = Kabs, KG = KG)
}

recursion_moments <- function(n, d) {
  K <- kernels(n, d)
  pn <- (2 / n)^d
  f  <- numeric(n); f[n] <- 1 - pn                       # mass in general states, by size
  J  <- list(); J[[1]] <- { z <- numeric(n); z[n] <- pn; z }   # tau = 0
  h <- 0
  while (sum(f) > 1e-16 && h < n) {
    J[[h + 2]] <- as.vector(f %*% K$Kabs)
    f <- as.vector(f %*% K$KG)
    h <- h + 1
  }
  J <- do.call(rbind, J)                                 # rows: tau = 0..h, cols: entry size m
  tau <- 0:(nrow(J) - 1)
  p_tau <- rowSums(J)
  E_tau <- sum(tau * p_tau); Var_tau <- sum(tau^2 * p_tau) - E_tau^2
  w_m   <- colSums(J)                                    # P(M_tau = m)
  E_R   <- sum(w_m * h_m)
  E_v   <- sum(w_m * v_m)                                # E{v(M_tau)}
  Var_h <- sum(w_m * h_m^2) - E_R^2                      # Var{h(M_tau)}
  Var_R <- E_v + Var_h
  Cov   <- sum(J * outer(tau, h_m)) - E_tau * E_R        # E(tau h(M_tau)) - E(tau) E(R)
  c(E_tau = E_tau, E_R = E_R, Var_tau = Var_tau, Var_R = Var_R, E_v = E_v, Var_h = Var_h, Cov = Cov,
    check_EH = E_tau + E_R - E_H, check_VH = Var_tau + Var_R + 2 * Cov - Var_H)
}

cat(sprintf("n = %d, E(H) = %.6f, Var(H) = %.6f\n", n, E_H, Var_H))
res <- NULL
for (d in d_set) {
  t0 <- Sys.time(); r <- recursion_moments(n, d)
  # dimension quantity at finite n (Definition 2; limiting allocation of Corollary 2)
  delta <- r[["E_tau"]] - (E_H - 2 / 3)
  D_tau <- r[["Var_tau"]] - (Var_H - 2 / 9)
  D_v   <- r[["E_v"]]
  D_h   <- r[["Var_h"]] - 2 / 9
  D_c   <- 2 * r[["Cov"]]
  fin   <- c(delta = delta, D_tau = D_tau, D_v = D_v, D_h = D_h, D_c = D_c, Var_R = r[["Var_R"]])
  lim   <- tab6[d, ]
  diff  <- fin - lim
  res <- rbind(res, data.frame(d = d, t(round(fin, 6)), max_abs_diff = max(abs(diff)),
                               diff_over_bound = max(abs(diff)) / (n^(-d) * log(n)),
                               sum_D = sum(fin[c("D_tau", "D_v", "D_h", "D_c")]),
                               check_EH = r[["check_EH"]], check_VH = r[["check_VH"]],
                               secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)))
  cat(sprintf("d = %2d  max|finite - limit| = %.2e  (bound ratio %.2f)  sum of D = %.1e  [%s s]\n",
              d, max(abs(diff)), max(abs(diff)) / (n^(-d) * log(n)), sum(fin[2:5]), res$secs[nrow(res)]))
}
print(res, digits = 7)
write.csv(res, file.path(OUT_DIR, sprintf("SuppTable6_zeta_limits_n%d.csv", n)), row.names = FALSE)
cat("\nExpected: d = 1 differs from the limit by O(log(n)/n) (ratio ~ 11);",
    "d >= 2 agree to about 6 decimals (Table 6 is rounded to 6 decimals).\n",
    "sum_D must be 0 and check_EH, check_VH must be ~1e-13 for every d.\n")
