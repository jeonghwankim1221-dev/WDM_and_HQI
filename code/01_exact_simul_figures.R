###############################################################################
# Figs. 2-3 and Extended Data Figs. 1-3
#
#  Fig. 2   allocation of the depth with dimension (exact curves + simulation)
#  Fig. 3   joint distribution of the entry depth and the depth, d = 1 and d -> Inf
#  Extended Data Fig. 2   composition of the observations by depth, d = 1 and d -> Inf
#  (Fig. 4, real data, is drawn by 05_table2_extdata_table2.R from the real-data run)
#  Extended Data Fig. 1   distribution of the entry size
#  Extended Data Fig. 3   depth-wise contributions to Var(tau) and Var(R)
#
#  Exact values : recursion of Corollary 1 (Theorem 4 kernels), d = Inf uses the
#                 limiting kernel of Corollary 2; residual law = Theorem 1.
#  Simulation   : basic model of Methods (iid Uniform, one variable per split,
#                 cut uniform between min and max), cached in SIM_RDS.
#
#  Output: OUT_DIR/Fig2, Fig3, ExtData_Fig1..ExtData_Fig3 (Fig. 1 is a schematic drawn separately) (.pdf/.png) and source data (source_data/*.csv)
###############################################################################

rm(list = ls())
suppressWarnings(suppressMessages(library(parallel)))

source(file.path(Sys.getenv("CODE_DIR", "."), "paths.R"))    # DATA_DIR (read-only), OUT_DIR (results)
SIM_RDS_IN  <- file.path(DATA_DIR, "sim_ch6.rds")            # cache shipped with the capsule
SIM_RDS_OUT <- file.path(OUT_DIR,  "sim_ch6.rds")            # written only when the simulation is rerun

CFG <- list(
  n_set    = c(256L, 512L),
  d_exact  = sort(unique(c(1:10, round(10^seq(1, 2, length.out = 12))))),
  d_sim    = c(1:10, 20L, 50L, 100L),
  n_trees  = 500L,          # trees per replicate
  n_rep    = 20L,           # replicates (a new sample each time)
  seed     = 20260928L,
  cores    = max(1L, floor(0.9 * parallel::detectCores())),
  run_sim  = TRUE,          # FALSE: exact curves only
  force_sim = FULL_RUN      # TRUE: ignore the cache and simulate again (set FULL=1)
)

# 1. exact computation #
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

depth_curves <- function(P) {                   # by depth h: P(tau>h), P(tau<=h<H), P(H>h), v_tau(h), v_R(h)
  Tm <- nrow(P); Rm <- ncol(P); L <- Tm + Rm
  t <- 0:(Tm - 1); r <- 0:(Rm - 1)
  Et <- sum(t * rowSums(P)); ER <- sum(r * colSums(P))
  gen <- ext <- alive <- vt <- vR <- numeric(L)
  for (h in 0:(L - 1)) {
    It <- outer(t > h, rep(TRUE, Rm), "&")              # tau > h (logical)
    IR <- outer(t, r, function(a, b) a <= h & h < a + b) # tau <= h < H
    gen[h + 1] <- sum(P[It]); ext[h + 1] <- sum(P[IR])
    alive[h + 1] <- sum(P[outer(t, r, function(a, b) h < a + b)])
    vt[h + 1] <- sum((outer(t, rep(1, Rm)) * P)[It]) - gen[h + 1] * Et
    vR[h + 1] <- sum((outer(rep(1, Tm), r) * P)[IR]) - ext[h + 1] * ER
  }
  data.frame(h = 0:(L - 1), gen = gen, ext = ext, alive = alive, v_tau = vt, v_R = vR)
}

heat_C <- function(J, ext, Kmax) {              # C[h+1, k+1] = P(tau = h, H = k)
  C <- matrix(0, nrow(J), Kmax + 1)
  for (m in seq_len(ncol(J))) {
    w <- J[, m]; if (!any(w > 0)) next
    e <- ext[[m]]
    for (h in which(w > 0) - 1) {
      k <- h + seq_along(e) - 1; ok <- k <= Kmax
      C[h + 1, k[ok] + 1] <- C[h + 1, k[ok] + 1] + w[h + 1] * e[ok]
    }
  }
  C
}

# ==================== 2. simulation ====================
one_tree <- function(Z, d) {                    # tau, M_tau, H for every observation of tree
  n <- nrow(Z)
  H <- integer(n); TAU <- rep(NA_integer_, n); MT <- rep(NA_integer_, n)
  stack <- list(list(idx = seq_len(n), h = 0L)); top <- 1L
  while (top > 0L) {
    nd <- stack[[top]]; stack[[top]] <- NULL; top <- top - 1L
    idx <- nd$idx; h <- nd$h; m <- length(idx)
    if (m == 1L) { H[idx] <- h; if (is.na(TAU[idx])) { TAU[idx] <- h; MT[idx] <- 1L }; next }
    S <- Z[idx, , drop = FALSE]
    v1 <- S[, 1L]; cand <- idx[(v1 == min(v1) | v1 == max(v1))]; cand <- cand[is.na(TAU[cand])]
    for (cc in cand) {
      row <- matrix(Z[cc, ], m, d, byrow = TRUE)
      if (all(colSums(S < row) == 0L | colSums(S > row) == 0L)) { TAU[cc] <- h; MT[cc] <- m }
    }
    repeat {
      q <- sample.int(d, 1L); v <- S[, q]; rg <- range(v)
      if (rg[2] > rg[1]) break
    }
    cut <- runif(1L, rg[1], rg[2]); L <- idx[v < cut]; R <- idx[v >= cut]
    top <- top + 1L; stack[[top]] <- list(idx = L, h = h + 1L)
    top <- top + 1L; stack[[top]] <- list(idx = R, h = h + 1L)
  }
  list(tau = TAU, M = MT, H = H)
}

one_rep <- function(seed, n, d, n_trees, keep_joint) {
  set.seed(seed)
  Z <- matrix(runif(n * d), n, d)
  tau <- M <- H <- integer(0)
  for (b in seq_len(n_trees)) { r <- one_tree(Z, d); tau <- c(tau, r$tau); M <- c(M, r$M); H <- c(H, r$H) }
  R <- H - tau
  mu <- function(x) mean(x); vr <- function(x) mean(x^2) - mean(x)^2
  gm <- tapply(R, M, mean); gv <- tapply(R, M, vr); gw <- tapply(R, M, length) / length(R)
  out <- c(E_tau = mu(tau), E_R = mu(R), E_H = mu(H), Var_tau = vr(tau), Var_R = vr(R),
           Cov = mean(tau * R) - mu(tau) * mu(R), Var_H = vr(H), E_M = mu(M),
           E_v = sum(gw * gv), Var_h = sum(gw * gm^2) - sum(gw * gm)^2)
  jt <- NULL
  if (keep_joint) jt <- table(factor(tau, levels = 0:80), factor(R, levels = 0:40))
  list(stats = out, joint = jt)
}

run_setting <- function(n, d, keep_joint) {
  seeds <- CFG$seed + 1000L * d + n + seq_len(CFG$n_rep)
  f <- function(s) one_rep(s, n, d, CFG$n_trees, keep_joint)
  if (.Platform$OS.type == "windows") {
    cl <- makeCluster(min(CFG$cores, CFG$n_rep))
    clusterExport(cl, c("one_tree", "one_rep", "CFG"), envir = globalenv())
    res <- parLapply(cl, seeds, function(s, n, d, nt, kj) one_rep(s, n, d, nt, kj),
                     n = n, d = d, nt = CFG$n_trees, kj = keep_joint)
    stopCluster(cl)
  } else res <- mclapply(seeds, f, mc.cores = min(CFG$cores, CFG$n_rep))
  S <- do.call(rbind, lapply(res, `[[`, "stats"))
  jt <- if (keep_joint) Reduce(`+`, lapply(res, `[[`, "joint")) else NULL
  list(mean = colMeans(S), se = apply(S, 2, sd) / sqrt(nrow(S)), joint = jt)
}

## ==================== run ====================
ext_by_n <- lapply(CFG$n_set, ext_dists); names(ext_by_n) <- CFG$n_set
cat("exact computation ...\n")
EX <- do.call(rbind, lapply(CFG$n_set, function(n)
  do.call(rbind, lapply(c(CFG$d_exact, Inf), function(d) exact_summary(n, d, ext_by_n[[as.character(n)]])))))
EX <- as.data.frame(EX)
write.csv(EX, file.path(OUT_DIR, "ExtDataTable1_exact_moments.csv"), row.names = FALSE)

SIM <- NULL
if (CFG$run_sim) {
  if (file.exists(SIM_RDS_IN) && !CFG$force_sim) { SIM <- readRDS(SIM_RDS_IN); cat("simulation loaded from cache\n") } else {
    SIM <- list()
    for (n in CFG$n_set) for (d in CFG$d_sim) {
      cat(sprintf("simulating n=%d d=%d\n", n, d))
      SIM[[paste(n, d)]] <- run_setting(n, d, keep_joint = (n == 256L && d %in% c(1L, 100L)))
    }
    saveRDS(SIM, SIM_RDS_OUT)
  }
}
sim_tab <- function(n, key) if (is.null(SIM)) NULL else
  sapply(CFG$d_sim, function(d) SIM[[paste(n, d)]]$mean[[key]])

## check against Extended Data Table 1 (n = 256)
chk <- subset(EX, n == 256 & (d %in% c(1:5, 10, 100) | is.infinite(d)))
print(round(chk[, c("d", "E_tau", "E_R", "E_M", "Var_tau", "Var_R", "Cov", "E_v", "Var_h")], 4))


## ==================== Config ====================
MM <- 1 / 25.4
W2 <- 180 * MM; W1 <- 88 * MM
PS <- 8
use_cairo <- isTRUE(capabilities("cairo"))
dev_both <- function(name, w, h, expr) {
  if (use_cairo) cairo_pdf(file.path(OUT_DIR, paste0(name, ".pdf")), width = w, height = h, pointsize = PS, family = "sans")
  else pdf(file.path(OUT_DIR, paste0(name, ".pdf")), width = w, height = h, pointsize = PS, family = "Helvetica")
  expr(); dev.off()
  png(file.path(OUT_DIR, paste0(name, ".png")), width = w, height = h, units = "in", res = 600, pointsize = PS,
      family = "sans", type = if (use_cairo) "cairo" else getOption("bitmapType"))
  expr(); dev.off()
}
SRC_DIR <- file.path(OUT_DIR, "source_data"); dir.create(SRC_DIR, showWarnings = FALSE, recursive = TRUE)
src_csv <- function(df, name) write.csv(df, file.path(SRC_DIR, paste0(name, ".csv")), row.names = FALSE)

col_tau <- "#0072B2"; col_R <- "#D55E00"; col_H <- "black"          # Okabe-Ito (colour-blind safe)
col_v <- "#0072B2"; col_h <- "#D55E00"; col_sim <- "#D55E00"; col_lim <- "#D55E00"
PT_fill <- 0.8; PT_open <- 1.1
pal_d <- c("#000000", "#0072B2", "#009E73", "#E69F00", "#CC79A7")   # d = 1, 2, 5, 10, 100
## panel label: bold lower-case letter; optional short title in plain text after it
panel_lab <- function(s, title = NULL) {
  lab <- if (is.null(title)) bquote(bold(.(s))) else bquote(bold(.(s)) ~~~~ .(title[[1]]))
  mtext(lab, side = 3, adj = 0, line = 0.4, cex = 1.1)
}
exn <- function(nn) EX[EX$n == nn & is.finite(EX$d), ]
xlab_d <- expression("Dimension " * italic(d) * " (log scale)")
xlab_h <- expression("Depth " * italic(h))
PAR0 <- list(mgp = c(2.1, 0.6, 0), tcl = -0.3, cex.axis = 0.9, cex.lab = 1, las = 1)

# ============================== Fig. 2 (fig1) ==============================
fig1 <- function() {
  do.call(par, PAR0)
  par(mfrow = c(2, 2), mar = c(3.6, 4.0, 1.6, 0.8), oma = c(2.2, 0, 0, 0))
  par(cex = 1)                                  # mfrow(2,2) would reduce text to 83%
  ltys <- c(1, 2); pchs <- c(16, 1)
  add_sim <- function(key, col) if (!is.null(SIM)) for (i in seq_along(CFG$n_set))
    points(CFG$d_sim, sim_tab(CFG$n_set[i], key), col = col, pch = pchs[i], cex = c(PT_fill, PT_open)[i], lwd = 1)
  EH_max <- max(EX$E_H[is.finite(EX$d)])
  # a
  plot(1, 1, log = "x", type = "n", xlim = c(1, 100), ylim = c(0, EH_max + 0.8), xlab = xlab_d, ylab = "Expected depth")
  for (i in seq_along(CFG$n_set)) { e <- exn(CFG$n_set[i])
    lines(e$d, e$E_tau, col = col_tau, lty = ltys[i], lwd = 1.2); lines(e$d, e$E_R, col = col_R, lty = ltys[i], lwd = 1.2)
    lines(e$d, e$E_H, col = col_H, lty = ltys[i], lwd = 1.2) }
  add_sim("E_tau", col_tau); add_sim("E_R", col_R); add_sim("E_H", col_H)
  legend("right", c(expression(E(tau)), expression(E(italic(R))), expression(E(italic(H)))),
         col = c(col_tau, col_R, col_H), lwd = 1.4, bty = "n", cex = 0.9)
  panel_lab("a")
  # b
  plot(1, 1, log = "x", type = "n", xlim = c(1, 100), ylim = c(0, 2.1), xlab = xlab_d, ylab = "Variance of residual depth")
  for (i in seq_along(CFG$n_set)) { e <- exn(CFG$n_set[i])
    lines(e$d, e$Var_R, lty = ltys[i], lwd = 1.2); lines(e$d, e$E_v, col = col_v, lty = ltys[i], lwd = 1.2)
    lines(e$d, e$Var_h, col = col_h, lty = ltys[i], lwd = 1.2) }
  add_sim("Var_R", "black"); add_sim("E_v", col_v); add_sim("Var_h", col_h)
  abline(h = 2 / 9, lty = 3)
  legend("topright", c(expression(Var(italic(R))), expression(E*"{"*italic(v)(italic(M)[tau])*"}"),
                       expression(Var*"{"*italic(h)(italic(M)[tau])*"}"), expression(italic(d) ~ "\u2192" ~ infinity * ":  2/9")),
         col = c("black", col_v, col_h, "black"), lty = c(1, 1, 1, 3), lwd = c(1.4, 1.4, 1.4, 1), bty = "n", cex = 0.9)
  panel_lab("b")
  # c
  plot(1, 1, log = "xy", type = "n", xlim = c(1, 100), ylim = c(1.5, 30), xlab = xlab_d,
       ylab = expression(E(italic(M)[tau]) * " (log scale)"))
  for (i in seq_along(CFG$n_set)) { e <- exn(CFG$n_set[i]); lines(e$d, e$E_M, lty = ltys[i], lwd = 1.2) }
  add_sim("E_M", "black"); abline(h = 5 / 3, lty = 3)
  legend("topright", expression(italic(d) ~ "\u2192" ~ infinity * ":  5/3"), lty = 3, bty = "n", cex = 0.9)
  panel_lab("c")
  # d
  plot(1, 1, log = "x", type = "n", xlim = c(1, 100), ylim = c(-2.1, 0.1), xlab = xlab_d, ylab = expression(Cov(tau * "," ~ italic(R))))
  for (i in seq_along(CFG$n_set)) { e <- exn(CFG$n_set[i]); lines(e$d, e$Cov, lty = ltys[i], lwd = 1.2) }
  add_sim("Cov", "black"); abline(h = 0, lty = 3)
  panel_lab("d")
  # common key
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE); plot.new()
  legend("bottom", c("Theory, n = 256", "Theory, n = 512", "Simulation, n = 256", "Simulation, n = 512"),
         lty = c(1, 2, NA, NA), pch = c(NA, NA, 16, 1), lwd = c(1.2, 1.2, NA, 1), pt.cex = c(1, 1, PT_fill, PT_open),
         horiz = TRUE, bty = "n", cex = 0.95, inset = 0.005)
}
dev_both("Fig2", W2, 150 * MM, fig1)
## source data: exact curves and simulated points
src_csv(EX[is.finite(EX$d), c("n", "d", "E_tau", "E_R", "E_H", "Var_R", "E_v", "Var_h", "E_M", "Cov")], "Fig2_exact")
if (!is.null(SIM)) {
  keys <- c("E_tau", "E_R", "E_H", "Var_R", "E_v", "Var_h", "E_M", "Cov")
  S1 <- do.call(rbind, lapply(CFG$n_set, function(nn) data.frame(n = nn, d = CFG$d_sim,
          sapply(keys, function(k) sim_tab(nn, k)))))
  src_csv(S1, "Fig2_simulation")
}

## ============================ Fig. 3 (fig2) ================================
n0 <- 256L; ext0 <- ext_by_n[["256"]]
J1 <- joint_tau_M(n0, 1); Jinf <- joint_tau_M(n0, Inf)
Kmax <- 25
C1 <- heat_C(J1, ext0, Kmax); Cinf <- heat_C(Jinf, ext0, Kmax)
PH <- bern_sum(2 / (2:n0))                          # Theorem 3, P(H = k), k = 0..n-1
EH_th <- sum((0:(n0 - 1)) * PH); VH_th <- sum((0:(n0 - 1))^2 * PH) - EH_th^2   # = 2H_n - 2, 2H_n + 2 - 4H_n^(2)
cat("column sums vs Theorem 3 (max abs diff): d=1", max(abs(colSums(C1) - PH[1:(Kmax + 1)])),
    " d=Inf", max(abs(colSums(Cinf) - PH[1:(Kmax + 1)])), "\n")
heat_cols <- hcl.colors(60, "YlOrRd", rev = TRUE)
zl <- c(-6, 0)
sim_joint <- function(dd) if (is.null(SIM)) NULL else SIM[[paste(256, dd)]]$joint

heat_block <- function(C, lab, title, sj) {
  Hh <- nrow(C) - 1; k <- 0:Kmax; h <- 0:Hh; ytop <- min(Hh, Kmax) + 0.5
  z <- log10(t(C)); z[!is.finite(z) | z < zl[1]] <- NA
  ph_s <- pt_s <- NULL
  if (!is.null(sj)) {
    tot <- outer(0:(nrow(sj) - 1), 0:(ncol(sj) - 1), `+`)
    ph_s <- tapply(as.vector(sj), as.vector(tot), sum); ph_s <- ph_s / sum(ph_s)
    pt_s <- rowSums(sj) / sum(sj)
  }
  # top strip: P(H = k)
  par(mar = c(0.2, 3.4, 1.6, 0.2))
  cs <- colSums(C)
  plot(k, cs, type = "n", xlim = c(-0.5, Kmax + 0.5), ylim = c(0, max(PH) * 1.15), axes = FALSE, xlab = "", ylab = "")
  rect(k - 0.4, 0, k + 0.4, cs, col = "grey75", border = NA)
  if (!is.null(ph_s)) { kk <- as.integer(names(ph_s)); okk <- kk <= Kmax
    points(kk[okk], ph_s[okk], pch = 16, cex = 0.85, col = col_sim) }
  text(Kmax + 0.5, max(PH) * 1.13, sprintf("mean %.3f\nsd %.3f", EH_th, sqrt(VH_th)), adj = c(1, 1), cex = 0.95)
  panel_lab(lab, title)
  # heat map
  par(mar = c(3.4, 3.4, 0.2, 0.2))
  image(k, h, z, zlim = zl, col = heat_cols, xlab = expression("Depth " * italic(k)),
        ylab = expression("Entry depth " * italic(h)), xlim = c(-0.5, Kmax + 0.5), ylim = c(-0.5, ytop))
  abline(0, 1, lty = 3)
  text(ytop - 1.4, ytop - 1.2, expression(paste(italic(k) == italic(h), "  (", italic(R) == 0, ")")), adj = c(1, 0.5), cex = 0.95)
  x0 <- 1; x1 <- 12; y0 <- ytop - 4.4; y1 <- ytop - 3.2
  xs <- seq(x0, x1, length.out = length(heat_cols) + 1)
  rect(xs[-length(xs)], y0, xs[-1], y1, col = heat_cols, border = NA)
  tk <- seq(zl[1], zl[2], by = 2); xt <- x0 + (tk - zl[1]) / diff(zl) * (x1 - x0)
  text(xt, y0 - 0.9, parse(text = paste0("10^", tk)), cex = 0.9)
  text(x0, y1 + 1.0, expression(P(tau == italic(h) * "," ~ italic(H) == italic(k))), adj = 0, cex = 0.95)
  box()
  # right strip: P(tau = h)
  par(mar = c(3.4, 0.2, 0.2, 0.4))
  rs <- rowSums(C)
  plot(rs, h, type = "n", ylim = c(-0.5, ytop), xlim = c(0, max(rs) * 1.12), axes = FALSE, xlab = "", ylab = "")
  rect(0, h - 0.4, rs, h + 0.4, col = "grey75", border = NA)
  if (!is.null(pt_s)) { hh <- 0:(length(pt_s) - 1); okh <- hh <= ytop
    points(pt_s[okh], hh[okh], pch = 16, cex = 0.85, col = col_sim) }
}
fig2 <- function() {
  do.call(par, PAR0)
  ## cells: a = 1 (top), 2 (map), 3 (right); b = 4, 5, 6; keys 7, 8
  layout(rbind(c(1, 7, 0, 4, 8), c(2, 3, 0, 5, 6)), widths = c(1, 0.24, 0.1, 1, 0.24), heights = c(0.3, 1))
  par(cex = 1)                                  # layout with two rows would reduce text
  heat_block(C1, "a", expression(italic(d) == 1), sim_joint(1))
  heat_block(Cinf, "b", expression(italic(d) ~ "\u2192" ~ infinity), sim_joint(100))
  par(mar = c(0, 0, 0, 0))
  for (cell in 1:2) { plot.new()
    legend("center", c("Exact", "Simulation"), fill = c("grey75", NA), border = NA,
           pch = c(NA, 16), col = c(NA, col_sim), pt.cex = 0.9, bty = "n", cex = 0.85, y.intersp = 1.3) }
}
dev_both("Fig3", W2, 95 * MM, fig2)
to_long <- function(C, dlab) { g <- expand.grid(h = 0:(nrow(C) - 1), k = 0:(ncol(C) - 1))
  data.frame(d = dlab, g, prob = as.vector(C)) }
src_csv(rbind(to_long(C1, "1"), to_long(Cinf, "Inf")), "Fig3_joint_tau_H")

# ============================== Extended Data Fig. 2 (fig3) ==============================
DC <- function(d) depth_curves(joint_tau_R(joint_tau_M(n0, d), ext0))
D1 <- DC(1); Dinf <- DC(Inf)
cat(sprintf("band areas: d=1 %.3f / %.3f, d=Inf %.3f / %.3f\n", sum(D1$gen), sum(D1$ext), sum(Dinf$gen), sum(Dinf$ext)))
hmax <- 20
fig3 <- function() {
  do.call(par, PAR0)
  par(mfrow = c(1, 2), mar = c(3.4, 3.8, 1.6, 0.8))
  band <- function(D, lab, title) {
    D <- D[D$h <= hmax, ]; h <- D$h
    plot(h, D$gen, type = "n", ylim = c(0, 1.08), yaxs = "i", xlab = xlab_h, ylab = "Proportion of observations")
    polygon(c(h, rev(h)), c(D$gen, rep(0, length(h))), col = "grey75", border = NA)
    polygon(c(h, rev(h)), c(D$gen + D$ext, rev(D$gen)), col = col_R, border = NA)
    lines(h, D$gen + D$ext, lwd = 0.7)
    box()
    legend("topright", c(sprintf("General state (area %.2f)", sum(D$gen)),
                         sprintf("Extreme state (area %.2f)", sum(D$ext))),
           fill = c("grey75", col_R), border = NA, bty = "n", cex = 0.85, inset = c(0.01, 0.01))
    panel_lab(lab, title)
  }
  band(D1, "a", expression(italic(d) == 1))
  band(Dinf, "b", expression(italic(d) ~ "\u2192" ~ infinity))
}
dev_both("ExtData_Fig2", W2, 72 * MM, fig3)
src_csv(rbind(data.frame(d = "1", D1[, c("h", "gen", "ext")]), data.frame(d = "Inf", Dinf[, c("h", "gen", "ext")])), "ExtData_Fig2_bands")

## ---- height limit (Discussion): quantities read from the bands of Extended Data Fig. 2 at depth 8 ----
##  P(H <= 8), P(tau > 8), and the areas of the two bands beyond depth 8, i.e. the expected
##  number of steps beyond the limit before entry, sum_{h>=8} P(tau > h), and after entry,
##  sum_{h>=8} P(tau <= h < H); the last column is the latter as a share of E(R) = sum_h P(tau <= h < H).
HL <- 8L
height_limit <- do.call(rbind, lapply(c(1, 2, 5, 10, Inf), function(dd) {
  D <- DC(dd); beyond <- D$h >= HL
  data.frame(d = as.character(dd), limit = HL,
             P_H_le_limit = 1 - D$alive[D$h == HL], P_tau_gt_limit = D$gen[D$h == HL],
             steps_beyond_before_entry = sum(D$gen[beyond]), steps_beyond_residual = sum(D$ext[beyond]),
             E_R = sum(D$ext))
}))
height_limit$residual_share_of_E_R <- height_limit$steps_beyond_residual / height_limit$E_R
print(height_limit, digits = 3)
src_csv(height_limit, "Discussion_height_limit")

# ============================== Extended Data Fig. 1 ==============================
sfig1 <- function() {
  do.call(par, PAR0)
  par(mar = c(3.4, 4.2, 0.8, 0.6))
  ds <- c(1, 2, 5, 10, 100)
  plot(1, 1, log = "xy", type = "n", xlim = c(1, n0), ylim = c(1e-6, 1), xlab = expression("Entry size " * italic(m) * " (log scale)"),
       ylab = expression(P(italic(M)[tau] == italic(m)) * " (log scale)"), yaxt = "n")
  axis(2, at = 10^(-6:0), labels = parse(text = paste0("10^", -6:0)))
  for (i in seq_along(ds)) {
    pm <- colSums(joint_tau_M(n0, ds[i])); m <- seq_len(n0); ok <- pm > 1e-7
    mb <- 2:(n0 - 1); bd <- 2 * (2 / mb)^ds[i] / (mb + 1); okb <- bd > 1e-7
    lines(mb[okb], bd[okb], col = pal_d[i], lty = 3, lwd = 0.8)
    lines(m[ok], pm[ok], col = pal_d[i], lwd = 1.2)
    ends <- c(min(m[ok]), if (ok[n0]) n0)
    points(ends, pm[ends], col = pal_d[i], pch = 16, cex = PT_fill)
  }
  points(c(1, 2), c(1 / 3, 2 / 3), pch = 4, col = col_lim, cex = 1.2, lwd = 1.6)
  legend("topright", c(paste("d =", ds), expression(italic(d) ~ "\u2192" ~ infinity), "Upper bound"),
         col = c(pal_d, col_lim, "grey30"), lty = c(rep(1, 5), NA, 3), pch = c(rep(NA, 5), 4, NA),
         lwd = c(rep(1.2, 5), 1.6, 0.8), bty = "n", cex = 0.8, inset = c(0.01, 0.01))
}
dev_both("ExtData_Fig1", W1, 75 * MM, sfig1)
src_csv(do.call(rbind, lapply(c(1, 2, 5, 10, 100), function(dd)
  data.frame(d = dd, m = seq_len(n0), prob = colSums(joint_tau_M(n0, dd))))), "ExtData_Fig1_entry_size")

# ============================== Extended Data Fig. 3 ==============================
sim_depth <- function(jt) { if (is.null(jt)) return(NULL); depth_curves(unclass(jt / sum(jt))) }
SD <- if (is.null(SIM)) list() else lapply(c(1, 100), function(dd) sim_depth(sim_joint(dd)))
sfig2 <- function() {
  do.call(par, PAR0)
  par(mar = c(3.4, 3.8, 0.8, 0.6))
  yr <- range(c(D1$v_tau, D1$v_R, Dinf$v_tau, Dinf$v_R)) * c(1, 1.15)
  plot(0, 0, type = "n", xlim = c(0, hmax), ylim = yr, xlab = xlab_h, ylab = "Contribution to variance")
  lines(D1$h, D1$v_tau, col = col_tau, lwd = 1.2); lines(Dinf$h, Dinf$v_tau, col = col_tau, lty = 2, lwd = 1.2)
  lines(D1$h, D1$v_R, col = col_R, lwd = 1.2); lines(Dinf$h, Dinf$v_R, col = col_R, lty = 2, lwd = 1.2)
  for (j in seq_along(SD)) { S <- SD[[j]]; if (is.null(S)) next
    ok <- S$h <= hmax & S$alive > 1e-3
    points(S$h[ok], S$v_tau[ok], col = col_tau, pch = c(16, 1)[j], cex = PT_fill)
    points(S$h[ok], S$v_R[ok], col = col_R, pch = c(16, 1)[j], cex = PT_fill) }
  abline(h = 0, lty = 3)
  legend("topright", c(expression(italic(v)[tau](italic(h))), expression(italic(v)[italic(R)](italic(h))),
                       expression(italic(d) == 1), expression(italic(d) ~ "\u2192" ~ infinity),
                       "Simulation, d = 1", "Simulation, d = 100"),
         col = c(col_tau, col_R, "black", "black", "black", "black"), lty = c(1, 1, 1, 2, NA, NA),
         pch = c(NA, NA, NA, NA, 16, 1), lwd = c(1.2, 1.2, 1.2, 1.2, NA, NA), bty = "n", cex = 0.8)
}
dev_both("ExtData_Fig3", W1, 75 * MM, sfig2)
src_csv(rbind(data.frame(d = "1", D1[, c("h", "v_tau", "v_R")]), data.frame(d = "Inf", Dinf[, c("h", "v_tau", "v_R")])), "ExtData_Fig3_variance_contributions")

cat("figures and source data written to", OUT_DIR, "\n")
