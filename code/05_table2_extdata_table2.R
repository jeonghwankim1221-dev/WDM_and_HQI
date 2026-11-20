###############################################################################
# Assembles Extended Data Table 2 and Table 2 from the *_components.csv files written
# by 04_realdata.R (rank-transformed data only is used for Extended Data Table 2's
# "rank-transformed" rows; raw data gives the "raw" rows; theory is common), and
# draws Fig. 4 from Extended Data Table 2.
# Under the default FULL=0 the real-data run uses reduced settings (8 subsamples x 100 trees),
# so the table and Fig. 4 differ from the paper (20 x 200) in the last digits; FULL=1 reproduces the paper.
###############################################################################

rm(list = ls())
.script_dir <- function() { a <- commandArgs(FALSE); k <- grep("^--file=", a); if (length(k)) return(dirname(normalizePath(sub("^--file=", "", a[k[1]])))); for (i in rev(seq_len(sys.nframe()))) { f <- sys.frame(i)$ofile; if (!is.null(f)) return(dirname(normalizePath(f))) }; getwd() }
source(file.path(Sys.getenv("CODE_DIR", .script_dir()), "paths.R"))

datasets <- c("quakes", "magic", "wine_white", "letter", "supercond")
read_one <- function(ds, mg) {
  f <- file.path(OUT_DIR, sprintf("realdata_%s_%s.csv", ds, mg))
  if (!file.exists(f)) stop("missing ", f, " -- run 04_realdata.R first")
  x <- read.csv(f, stringsAsFactors = FALSE); x
}

# ---------- Extended Data Table 2
## mean depth (tau, R) and variance of depth (tau, R), theory / raw / rank, with percent deviation from theory.
rows62 <- list()
for (ds in datasets) {
  mr <- read_one(ds, "rank"); mo <- read_one(ds, "raw")
  mmom_r <- mr[mr$kind == "moment", ]; mmom_o <- mo[mo$kind == "moment", ]
  gth <- function(m, nm) m$theory[m$name == nm][1]
  gem <- function(m, nm) m$empirical[m$name == nm][1]
  th  <- c(E_tau = gth(mmom_r, "E_tau"), E_R = gth(mmom_r, "E_R"),
           V_tau = gth(mmom_r, "Var_tau"), V_R = gth(mmom_r, "Vh") + gth(mmom_r, "Ev"))
  # Var(R) = Ev (within entry sizes) + Vh (between entry sizes); see Methods
  raw <- c(E_tau = gem(mmom_o, "E_tau"), E_R = gem(mmom_o, "E_R"),
           V_tau = gem(mmom_o, "Var_tau"), V_R = gem(mmom_o, "Ev") + gem(mmom_o, "Vh"))
  rnk <- c(E_tau = gem(mmom_r, "E_tau"), E_R = gem(mmom_r, "E_R"),
           V_tau = gem(mmom_r, "Var_tau"), V_R = gem(mmom_r, "Ev") + gem(mmom_r, "Vh"))
  pct <- function(v) 100 * (v - th) / th
  na4 <- setNames(rep(NA_real_, 4), paste0(names(th), "_pct"))
  rows62[[length(rows62) + 1]] <- data.frame(dataset = ds, kind = "theory",
                                              as.list(th), as.list(na4), check.names = FALSE)
  rows62[[length(rows62) + 1]] <- data.frame(dataset = ds, kind = "raw",
                                              as.list(raw), as.list(setNames(pct(raw), names(na4))), check.names = FALSE)
  rows62[[length(rows62) + 1]] <- data.frame(dataset = ds, kind = "rank_transformed",
                                              as.list(rnk), as.list(setNames(pct(rnk), names(na4))), check.names = FALSE)
}
TAB62 <- do.call(rbind, rows62)
write.csv(TAB62, file.path(OUT_DIR, "ExtDataTable2_real_data_depths.csv"), row.names = FALSE)

# ---------- Table 2
## dimension quota (theory D at finite d) and arrangement quota (empirical - theory D), for the four variance components and the mean (delta); 
## rank-transformed data, as used in Results (real data).
rows71 <- list()
for (ds in datasets) {
  x <- read_one(ds, "rank"); Dt <- x[x$kind == "D", ]
  dim_q <- setNames(Dt$theory, Dt$name); arr_q <- setNames(Dt$empirical - Dt$theory, Dt$name)
  rows71[[length(rows71) + 1]] <- data.frame(dataset = ds, quota = "dimension",   t(dim_q))
  rows71[[length(rows71) + 1]] <- data.frame(dataset = ds, quota = "arrangement", t(arr_q))
}
TAB71 <- do.call(rbind, rows71)
write.csv(TAB71, file.path(OUT_DIR, "Table2_dimension_arrangement.csv"), row.names = FALSE)

cat("Extended Data Table 2:\n"); print(TAB62, row.names = FALSE)
cat("\nTable 2:\n"); print(TAB71, row.names = FALSE)
cat("\nwritten:", file.path(OUT_DIR, "ExtDataTable2_real_data_depths.csv"), "and Table2_dimension_arrangement.csv\n")

# ---------- Fig. 4
## mean and variance of the depth split at the entry depth, theory vs rank-transformed data, drawn from Extended Data Table 2 above.
MM <- 1 / 25.4; W2 <- 180 * MM; PS <- 8
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
PAR0 <- list(mgp = c(2.1, 0.6, 0), tcl = -0.3, cex.axis = 0.9, cex.lab = 1, las = 1)
panel_lab <- function(s) mtext(bquote(bold(.(s))), side = 3, adj = 0, line = 0.4, cex = 1.1)

## reference lines: E(H) and Var(H) for the subsample size psi (Theorem 3), independent of d
psi <- 256L                                             # CFG$psi of 04_realdata.R
Hn  <- sum(1 / seq_len(psi)); Hn2 <- sum(1 / seq_len(psi)^2)
EH_ref <- 2 * Hn - 2; VH_ref <- 2 * Hn + 2 - 4 * Hn2

RD <- TAB62[TAB62$kind %in% c("theory", "rank_transformed"), c("dataset", "kind", "E_tau", "E_R", "V_tau", "V_R")]
RD$kind <- ifelse(RD$kind == "theory", "theory", "rank")
DNAME <- c(quakes = "Quakes", magic = "MAGIC", wine_white = "White Wine", letter = "Letter", supercond = "Supercond.")
DDIM  <- c(quakes = 5, magic = 10, wine_white = 11, letter = 16, supercond = 81)
fig4 <- function() {
  do.call(par, PAR0)
  par(mfrow = c(1, 2), mar = c(3.6, 3.8, 1.6, 0.6))
  c_lo <- c(theory = "grey75", rank = "#56B4E9")          # before entry (tau): light
  c_hi <- c(theory = "grey30", rank = "#0072B2")          # after entry (R): dark
  X <- RD; nk <- 2; gap <- 3
  pos <- as.vector(outer(0:(nk - 1), gap * (0:4), `+`)) + 1
  grp <- function(which_mean) {
    top <- if (which_mean) X$E_tau + X$E_R else X$V_tau + X$V_R
    low <- if (which_mean) X$E_tau else X$V_tau
    ref <- if (which_mean) EH_ref else VH_ref
    plot(0, 0, type = "n", xlim = c(0.3, max(pos) + 0.7), ylim = c(0, max(top) * 1.4), yaxs = "i", xaxt = "n",
         xlab = "", ylab = if (which_mean) "Mean depth" else "Variance of depth")
    rect(pos - 0.4, 0, pos + 0.4, low, col = c_lo[X$kind], border = NA)
    rect(pos - 0.4, low, pos + 0.4, top, col = c_hi[X$kind], border = NA)
    abline(h = ref, lty = 3)
    ctr <- gap * (0:4) + (nk + 1) / 2
    axis(1, at = ctr, labels = FALSE, tcl = -0.2)
    axis(1, at = ctr, labels = DNAME[unique(X$dataset)], tick = FALSE, line = -0.4, cex.axis = 0.85)
    axis(1, at = ctr, labels = paste0("(d = ", DDIM[unique(X$dataset)], ")"), tick = FALSE, line = 0.4, cex.axis = 0.8)
    q <- c("\u03c4", "R")
    labs <- c(paste0(q[1], ", theory"), paste0(q[2], ", theory"), paste0(q[1], ", rank-transformed data"),
              paste0(q[2], ", rank-transformed data"), sprintf("H, theory (%.3f)", ref))
    legend("topright", labs, fill = c(c_lo["theory"], c_hi["theory"], c_lo["rank"], c_hi["rank"], NA),
           border = NA, lty = c(NA, NA, NA, NA, 3), bty = "n", cex = 0.8, inset = c(0.01, 0.01))
    panel_lab(if (which_mean) "a" else "b")
  }
  grp(TRUE); grp(FALSE)
}
dev_both("Fig4", W2, 80 * MM, fig4)
write.csv(transform(RD, data = DNAME[dataset], d = DDIM[dataset]), file.path(SRC_DIR, "Fig4_real_data.csv"), row.names = FALSE)
cat(sprintf("Fig. 4 written from Extended Data Table 2 (reference E(H) = %.3f, Var(H) = %.3f for psi = %d)\n", EH_ref, VH_ref, psi))
cat("\nwritten:", file.path(OUT_DIR, "ExtDataTable2_real_data_depths.csv"), "and Table2_dimension_arrangement.csv\n")
