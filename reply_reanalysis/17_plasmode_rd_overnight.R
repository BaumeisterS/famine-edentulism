# =====================================================================
#  17_plasmode_rd_overnight.R — extended plasmode simulation (overnight)
#
#  Adds the RD estimators the first simulation (04/04c) left out, and new
#  data-generating processes for sensitivity:
#   RD estimators per replication
#    - conventional local-linear RDiT at fixed bandwidths 9.52 (published),
#      7, 5, 3.84 (the honest fit's bandwidth) and 3 years; Donut1
#    - rdrobust exactly as published (covariates, masspoints adjust) with the
#      MSE-optimal bandwidth RE-SELECTED in every replication: conventional,
#      bias-corrected and robust inference
#    - RDHonest (M by rule of thumb, re-estimated every replication), on the
#      raw pooled outcome (as published) and on the study-demeaned outcome
#    - plus the reduced form, local ITS and dose-response (for the GLS)
#   DGPs (null = the logistic fit of 04/04c, no famine term)
#    - base: independent outcomes around smooth study-specific cohort trends
#    - shock: base + study x birth-year random shocks on the logit scale with
#      the SD estimated from the data (cohort-level shocks that the base DGP
#      lacks and province/cell clustering may not capture)
#   effect shapes (mean effect among the exposed = delta unless noted)
#    U, G, R (as 04), S1, S2 (as 04c), G_RR (graded, proportional to the
#    baseline risk = multiplicative effect), G_nonCN (graded effect only in
#    the non-Chinese famines; delta = mean among the exposed in those famines)
#
#  Resumable: every block of replications is saved to output/17_chunks/;
#  re-running skips finished blocks. Stops launching blocks after
#  SIM_DEADLINE. Summary: 17b_plasmode_rd_summary.R (sourced at the end).
#  Seeds: 20261114 + 7919 x cell index + replication (deterministic).
# =====================================================================
source("01_common.R")
Sys.setenv(REANALYSIS_DIR = RS)   # inherited by the parallel workers
suppressPackageStartupMessages(library(parallel))
CH <- file.path(OUT, "17_chunks"); dir.create(CH, showWarnings = FALSE)
LOG <- file.path(OUT, "17_progress.log")
logf <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
DEADLINE <- as.POSIXct(Sys.getenv("SIM_DEADLINE", "2100-01-01 00:00:00"))
NW <- as.integer(Sys.getenv("SIM_WORKERS", 12)); BLOCK <- as.integer(Sys.getenv("SIM_BLOCK", 36))
logf("start: workers %d, block %d, deadline %s", NW, BLOCK, format(DEADLINE))

# ---------------- set-up (identical on master and workers) -----------
setup <- function() {
  source(file.path(Sys.getenv("REANALYSIS_DIR"), "01_common.R"))
  suppressPackageStartupMessages({ library(rdrobust); library(RDHonest) })
  d <<- make_analytic(load_pooled(4)); NN <<- nrow(d)
  S <- do.call(cbind, lapply(sort(unique(d$study)), function(k) { M <- matrix(0, NN, 6); r <- d$study == k; M[r, ] <- ns(d$byear[r], df = 6); M }))
  g0 <- glm.fit(drop_collinear(cbind(1, S, cov_block(d))), d$edent, family = binomial())
  p0 <<- as.numeric(g0$fitted.values); eta0 <<- as.numeric(g0$linear.predictors)
  cellf <<- as.integer(factor(paste(d$study, d$byear))); ncell <<- max(cellf)
  ex <- d$exposed15 == 1
  wdose <- function(wfun) { out <- numeric(NN); for (a in 0:15) out <- out + wfun(a) * overlap_years(d$byear, d$famine_start, d$famine_end, a, a); out }
  G <- d$dose15 / mean(d$dose15[ex]); nonCN <- !(d$study %in% 1:2)
  SH <- list(U = as.numeric(ex), G = G,
             R = { z <- d$dose15 + 0.5 * d$rec15; z / mean(z[ex]) },
             S1 = wdose(function(a) ifelse(a <= 5, 1, pmax(0, (12 - a) / 7))),
             S2 = wdose(function(a) as.numeric(a <= 11)),
             G_RR = { z <- G * p0; z / mean(z[ex]) },
             G_nonCN = { z <- G * nonCN; z / mean(z[ex & nonCN]) })
  SH$S1 <- SH$S1 / mean(SH$S1[ex]); SH$S2 <- SH$S2 / mean(SH$S2[ex])
  SHAPES <<- SH
  md <- mean(d$dose15[ex])
  CHEAP <<- list(RDiT = make_rd(d), RDiT_h7 = make_rd(d, h = 7), RDiT_h5 = make_rd(d, h = 5),
                 RDiT_h384 = make_rd(d, h = 3.84), RDiT_h3 = make_rd(d, h = 3), Donut1 = make_rd(d, donut = 1),
                 RF = make_rf(d), ITS_local = make_its(d, weighted = TRUE, win = 15, trend = "study", centre = "cutoff"),
                 Dose = make_dose(d, win = 15, df = 4))
  CSC <<- c(RDiT = 1, RDiT_h7 = 1, RDiT_h5 = 1, RDiT_h384 = 1, RDiT_h3 = 1, Donut1 = 1, RF = 1, ITS_local = 1, Dose = md)
  CM <<- model.matrix(~ male + factor(father_ses) + factor(study), data = d)[, -1]
  invisible(TRUE)
}
all_estimators <- function(y) {
  ch <- unlist(lapply(names(CHEAP), function(e) { r <- CHEAP[[e]](y) * CSC[e]; setNames(as.numeric(r), paste0(e, c("_est", "_se"))) }))
  rr <- if (Sys.getenv("SIM_FAST") == "1") NULL else tryCatch(suppressWarnings(rdrobust(y, d$run, c = 0, p = 1, kernel = "triangular", bwselect = "mserd", covs = CM, masspoints = "adjust")),
                 error = function(e) NULL)
  rd <- if (is.null(rr)) setNames(rep(NA_real_, 6), c("rdr_conv_est", "rdr_conv_se", "rdr_bc_est", "rdr_rob_se", "rdr_h", "rdr_b")) else
    c(rdr_conv_est = rr$coef[1, 1], rdr_conv_se = rr$se[1, 1], rdr_bc_est = rr$coef[2, 1], rdr_rob_se = rr$se[3, 1],
      rdr_h = rr$bws[1, 1], rdr_b = rr$bws[2, 1])
  hon <- function(yy, lab) {
    f <- tryCatch(suppressWarnings(suppressMessages(RDHonest(yy ~ run, data = data.frame(yy = yy, run = d$run), cutoff = 0,
                  kern = "triangular", opt.criterion = "MSE", sclass = "H"))), error = function(e) NULL)
    v <- if (is.null(f)) rep(NA_real_, 7) else { c1 <- as.data.frame(f$coefficients)
      c(c1$estimate[1], c1$std.error[1], c1$conf.low[1], c1$conf.high[1], c1$conf.high.onesided[1], c1$bandwidth[1], c1$M[1]) }
    setNames(v, paste0(lab, c("_est", "_se", "_lb", "_ub", "_ub1", "_h", "_M")))
  }
  c(ch, rd, hon(y, "hon_raw"), hon(y - ave(y, d$study), "hon_dm"))
}
one_rep <- function(task) {
  set.seed(task$seed)
  p <- if (task$dgp == "shock") plogis(eta0 + rnorm(ncell, 0, TAU)[cellf]) else p0
  dp <- if (task$delta > 0) task$delta * SHAPES[[task$shape]] else 0
  y <- rbinom(NN, 1, pmin(pmax(p + dp, 0), 1))
  c(rep = task$rep, all_estimators(y))
}

setup()
# cohort-shock SD (study x birth-year cells) beyond the smooth null fit
TAU <- tryCatch({
  suppressPackageStartupMessages(library(lme4))
  m <- glmer(edent ~ 1 + (1 | cellf), data = data.frame(edent = d$edent, cellf = factor(cellf)), offset = eta0, family = binomial(),
             control = glmerControl(optimizer = "bobyqa"))
  as.numeric(attr(VarCorr(m)$cellf, "stddev"))
}, error = function(e) { logf("glmer failed (%s); using SD 0.15", conditionMessage(e)); 0.15 })
TAU_EST <- TAU
TAU <- as.numeric(Sys.getenv("SIM_TAU", "0.15"))   # stress-test SD used in the reported run
logf("cohort-shock SD on the logit scale (study x birth-year): estimated %.3f (glmer; moment estimate also 0, Pearson dispersion 0.82 over 455 cells with n >= 50); used in the shock DGP: %.3f (stress test)", TAU_EST, TAU)
if (Sys.getenv("SIM_TEST") == "1") {      # dry run: one replication per scenario type, cheap estimators only
  Sys.setenv(SIM_FAST = "1")
  for (k in list(list("base", "U", 0), list("shock", "G", .02), list("base", "S1", .02), list("base", "S2", .02),
                 list("base", "G_RR", .02), list("base", "G_nonCN", .02), list("base", "R", .02))) {
    v <- one_rep(list(rep = 1, seed = 1, dgp = k[[1]], shape = k[[2]], delta = k[[3]]))
    cat(sprintf("%-6s %-8s RDiT %+.3f  RDiT_h384 %+.3f  RF %+.3f  ITS %+.3f  Dose %+.3f  (n values %d)\n", k[[1]], k[[2]],
                100 * v["RDiT_est"], 100 * v["RDiT_h384_est"], 100 * v["RF_est"], 100 * v["ITS_local_est"], 100 * v["Dose_est"], length(v)))
  }
  print(sapply(SHAPES, function(z) c(mean_exposed = mean(z[d$exposed15 == 1]), share_pos = mean(z > 0))))
  quit(save = "no")
}
obs <- all_estimators(d$edent)
true_mean <- lapply(SHAPES, function(z) mean(z[d$exposed15 == 1]))   # pooled mean effect among all exposed per unit delta
saveRDS(list(obs = obs, TAU = TAU, TAU_EST = TAU_EST, true_mean = true_mean), file.path(OUT, "17_observed.rds"))
logf("observed estimators computed")

# ---------------- scenarios and block schedule -----------------------
SC <- data.frame(cell = c("null", "U2", "G2", "R2", "S1_2", "S2_2", "shock_null", "shock_G2", "G_RR2", "G_nonCN2", "G1", "S2_3"),
                 dgp = c("base", "base", "base", "base", "base", "base", "shock", "shock", "base", "base", "base", "base"),
                 shape = c("U", "U", "G", "R", "S1", "S2", "U", "G", "G_RR", "G_nonCN", "G", "S2"),
                 delta = c(0, .02, .02, .02, .02, .02, 0, .02, .02, .02, .01, .03), stringsAsFactors = FALSE)
SC$id <- seq_len(nrow(SC))
# priority: every core cell to 144, then sensitivity cells to 144, then all to 288, then null to 576, then extras
targets <- list(list(cells = SC$cell[1:6], to = 144), list(cells = SC$cell[7:10], to = 144),
                list(cells = SC$cell[1:6], to = 288), list(cells = SC$cell[7:10], to = 288),
                list(cells = "null", to = 576), list(cells = c("shock_null"), to = 432),
                list(cells = SC$cell[11:12], to = 144), list(cells = SC$cell[2:6], to = 432))
# default 4 = the reported run (288 replications in each of the 10 cells); SIM_MAX_TARGET = 8 adds the optional cells
targets <- targets[seq_len(min(length(targets), as.integer(Sys.getenv("SIM_MAX_TARGET", "4"))))]
blocks <- do.call(rbind, lapply(targets, function(t) do.call(rbind, lapply(t$cells, function(cl) {
  st <- seq(1, t$to, by = BLOCK); data.frame(cell = cl, from = st, to = pmin(st + BLOCK - 1, t$to)) }))))
blocks <- blocks[!duplicated(blocks[, c("cell", "from")]), ]
fn <- function(b) file.path(CH, sprintf("%s_%04d-%04d.rds", b$cell, b$from, b$to))
logf("%d blocks scheduled (%d replications); %d already done", nrow(blocks), sum(blocks$to - blocks$from + 1),
     sum(file.exists(sapply(seq_len(nrow(blocks)), function(i) fn(blocks[i, ])))))

cl <- makeCluster(NW)
clusterExport(cl, c("setup", "all_estimators", "one_rep", "TAU"))
invisible(clusterEvalQ(cl, setup()))
logf("workers ready")
t_block <- NA
for (i in seq_len(nrow(blocks))) {
  b <- blocks[i, ]; f <- fn(b); if (file.exists(f)) next
  if (!is.na(t_block) && Sys.time() + t_block > DEADLINE) { logf("deadline reached; stopping before block %s %d-%d", b$cell, b$from, b$to); break }
  s <- SC[SC$cell == b$cell, ]
  tasks <- lapply(b$from:b$to, function(r) list(rep = r, seed = 20261114 + 7919 * s$id + r, dgp = s$dgp, shape = s$shape, delta = s$delta))
  t1 <- Sys.time()
  res <- parLapplyLB(cl, tasks, one_rep)
  saveRDS(do.call(rbind, res), f)
  t_block <- difftime(Sys.time(), t1, units = "secs")
  logf("block %-10s %4d-%4d done in %.1f min", b$cell, b$from, b$to, as.numeric(t_block) / 60)
}
stopCluster(cl)
logf("simulation loop finished")
source(file.path(RS, "17b_plasmode_rd_summary.R"))
logf("DONE 17")
