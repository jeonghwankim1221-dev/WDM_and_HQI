###############################################################################
# Checks that the closed forms of Theorem 5 and Appendix D agree with the
# recursion of Corollary 1 (the values in exact_ch6.csv, written by 01_exact_simul_figures.R).
#
#  c_n(m)  = 2n / {m(m+1)}   for m < n,   c_n(n) = 1             (Lemma 7)
#  s(m)    = 2^d / m^(d-1)   (expected number of extreme observations in a node of size m)
#  h(m)    = H_{m-1}                                             (Theorem 1)
#  D_n(m)  = 2(H_n - H_{m+1}) + 1  for m < n,   D_n(n) = 0        (Appendix C.2)
#  E(R)         = (1/n) sum_m c_n(m) s(m)                   (Theorem 5; Lemma 6)
#  E(R^2)       = (1/n) sum_m c_n(m) s(m) (2 h(m) - 1)
#  Var(R)       = E(R^2) - E(R)^2
#  E(H R)       = (1/n) sum_m c_n(m) s(m) {D_n(m) + h(m)}
#  Cov(tau, R)  = E(H R) - E(H) E(R) - Var(R)
#  Var(tau)     = Var(H) - Var(R) - 2 Cov(tau, R)
#  E(H) = 2 H_n - 2 and Var(H) = 2 H_n + 2 - 4 H_n^(2) are exact for every d (Theorem 3) and are not re-derived here.
#
# Output: OUT_DIR/closed_form_check.csv, and the maximum absolute difference
###############################################################################

rm(list = ls())
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))

EXF <- file.path(OUT_DIR, "ExtDataTable1_exact_moments.csv")
if (!file.exists(EXF)) stop("missing ", EXF, " -- run 01_exact_simul_figures.R first")
EX <- read.csv(EXF)
EX <- EX[is.finite(EX$d), ]                     # drop the d = Inf row(s); closed forms below are for finite d

Hn  <- function(n) sum(1 / seq_len(n))
Hn2 <- function(n) sum(1 / seq_len(n)^2)

closed_form <- function(n, d) {
  m  <- 2:n                                     # includes the root, m = n
  cn <- ifelse(m == n, 1, 2 * n / (m * (m + 1))) # Lemma 7; c_n(n) = 1 (the root itself)
  s  <- (2 / m)^d * m                           # s(m) = 2^d / m^(d-1)
  h  <- sapply(m, function(k) Hn(k - 1))        # h(m) = H_{m-1}
  Dn <- ifelse(m == n, 0, 2 * (Hn(n) - sapply(pmin(m + 1, n), Hn)) + 1)   # Appendix C.2; D_n(n) = 0 (root, depth 0)

  E_R   <- sum(cn * s) / n                    # Theorem 5: E(R) = (1/n) sum c_n(m) s(m)
  E_R2  <- sum(cn * s * (2 * h - 1)) / n
  Var_R <- E_R2 - E_R^2
  E_HR  <- sum(cn * s * (Dn + h)) / n

  E_H   <- 2 * Hn(n) - 2                        # Theorem 3, exact for every d
  Var_H <- 2 * Hn(n) + 2 - 4 * Hn2(n)           # Theorem 3, exact for every d

  Cov_tau_R <- E_HR - E_H * E_R - Var_R
  Var_tau   <- Var_H - Var_R - 2 * Cov_tau_R
  c(E_R = E_R, Var_R = Var_R, Cov = Cov_tau_R, Var_tau = Var_tau)
}

rows <- list()
for (i in seq_len(nrow(EX))) {
  n <- EX$n[i]; d <- EX$d[i]
  cf <- closed_form(n, d)
  rec <- c(E_R = EX$E_R[i], Var_R = EX$Var_R[i], Cov = EX$Cov[i], Var_tau = EX$Var_tau[i])
  rows[[i]] <- data.frame(n = n, d = d, quantity = names(cf),
                           closed_form = unname(cf), recursion = unname(rec),
                           abs_diff = abs(unname(cf) - unname(rec)))
}
TAB <- do.call(rbind, rows)
write.csv(TAB, file.path(OUT_DIR, "closed_form_check.csv"), row.names = FALSE)

cat(sprintf("compared %d (n, d, quantity) triples\n", nrow(TAB)))
cat(sprintf("maximum absolute difference between closed form and recursion: %.3e\n", max(TAB$abs_diff)))
print(aggregate(abs_diff ~ quantity, TAB, max))
