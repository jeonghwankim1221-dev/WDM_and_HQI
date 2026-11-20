###############################################################################
# 11_supp_note11_third_moment.R -- Supplementary Note 11, Proposition S4:
#   limit of E(R^3) as n -> Inf, closed form in single zeta values against an
#   independent numerical evaluation (base R, double precision).
#
#   closed form : mu3 = rho1 + 6 * 2^(d+1) * { sum_{j=2}^{d} (-1)^(d-j) zeta(j,1,1) + (-1)^(d-1) },
#                 rho1 = 2^(d+1) * { sum_{j=2}^{d} (-1)^(d-j) (zeta(j)-1) + (-1)^(d-1)/2 },
#                 zeta(j,1,1) reduced to single zeta values by the generating
#                 function of Borwein, Bradley and Broadhurst (Supplementary Eq. (35)).
#   numerical   : sum_{m>=2} w_m {1 + 6 e_1(m) + 6 e_2(m)},  w_m = 2^(d+1)/(m^d (m+1)),
#                 e_1(m) = H_{m-1} - 1,  e_2(m) = {(H_{m-1}-1)^2 - (H^(2)_{m-1} - 1)}/2,
#                 summed directly to m = K and continued by an Euler-Maclaurin tail
#                 with H_{x-1} = digamma(x) + gamma and H^(2)_{x-1} = zeta(2) - trigamma(x).
#   output      : OUT_DIR/third_moment_check.csv (d, closed form, numerical, difference)
###############################################################################
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))

## single zeta values from the polygamma function: zeta(k) = (-1)^k psi^(k-1)(1) / (k-1)!
zeta1 <- function(k) (-1)^k * psigamma(1, deriv = k - 1) / factorial(k - 1)

## bivariate power series, total degree <= N, stored as (N+1) x (N+1) matrices [x^a, y^b]
N <- 14
series_mult <- function(A, B) {
  C <- matrix(0, N + 1, N + 1)
  for (a in 0:N) for (b in 0:(N - a)) if (A[a + 1, b + 1] != 0)
    for (c in 0:(N - a - b)) for (e in 0:(N - a - b - c)) C[a + c + 1, b + e + 1] <- C[a + c + 1, b + e + 1] + A[a + 1, b + 1] * B[c + 1, e + 1]
  C
}
## S(x, y) = sum_{k>=2} zeta(k)/k * (x^k + y^k - (x+y)^k)
S <- matrix(0, N + 1, N + 1)
for (k in 2:N) for (a in 0:k) { b <- k - a
  S[a + 1, b + 1] <- S[a + 1, b + 1] + zeta1(k) / k * ((b == 0) + (a == 0) - choose(k, a)) }
## 1 - exp(S) = -(S + S^2/2! + S^3/3! + ...)
G <- matrix(0, N + 1, N + 1); P <- diag(0, N + 1); P[1, 1] <- 1
for (n in 1:N) { P <- series_mult(P, S); G <- G + P / factorial(n) }
G <- -G
zeta_j11 <- function(j) G[j - 1 + 1, 3 + 1]            # coefficient of x^(j-1) y^3  (Supplementary Eq. (35), a = j-2, b = 2)
stopifnot(abs(zeta_j11(2) - zeta1(4)) < 1e-12, abs(zeta_j11(3) - (2 * zeta1(5) - zeta1(2) * zeta1(3))) < 1e-12)

closed_mu3 <- function(d) {
  s    <- if (d >= 2) sum(sapply(2:d, function(j) (-1)^(d - j) * (zeta1(j) - 1))) else 0
  rho1 <- 2^(d + 1) * (s + (-1)^(d - 1) / 2)
  z    <- if (d >= 2) sum(sapply(2:d, function(j) (-1)^(d - j) * zeta_j11(j))) else 0
  rho1 + 6 * 2^(d + 1) * (z + (-1)^(d - 1))
}

numeric_mu3 <- function(d, K = 2e5) {
  m  <- 2:K
  H1 <- cumsum(1 / (1:K))[m - 1]; H2 <- cumsum(1 / (1:K)^2)[m - 1]     # H_{m-1}, H^(2)_{m-1}
  w  <- 2^(d + 1) / (m^d * (m + 1))
  f  <- function(x) { h1 <- digamma(x) - digamma(1); h2 <- zeta1(2) - trigamma(x)
    2^(d + 1) / (x^d * (x + 1)) * (1 + 6 * (h1 - 1) + 3 * ((h1 - 1)^2 - (h2 - 1))) }
  direct <- sum(w * (1 + 6 * (H1 - 1) + 3 * ((H1 - 1)^2 - (H2 - 1))))
  g  <- function(u) { v <- f(K * exp(u)) * K * exp(u); v[!is.finite(v)] <- 0; v }   # tail with x = K e^u: smooth, exponentially decaying
  I  <- integrate(g, 0, Inf, rel.tol = 1e-11, subdivisions = 1000L)$value
  h  <- K * 1e-4; f1 <- (f(K + h) - f(K - h)) / (2 * h)                  # Euler-Maclaurin: - f(K)/2 - B2/2! f'(K)
  direct + I - f(K) / 2 - f1 / 12
}

res <- do.call(rbind, lapply(1:8, function(d) { a <- closed_mu3(d); b <- numeric_mu3(d)
  data.frame(d = d, closed_form = a, numerical = b, abs_diff = abs(a - b), rel_diff = abs(a - b) / abs(a)) }))
print(format(res, digits = 15))
cat(sprintf("\nmax relative difference: %.2e\n", max(res$rel_diff)))
stopifnot(abs(res$closed_form[1] - 26) < 1e-12)                          # d = 1: E(R^3) -> 26
write.csv(res, file.path(OUT_DIR, "third_moment_check.csv"), row.names = FALSE)
