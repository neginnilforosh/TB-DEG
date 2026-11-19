## Soft-thresholding record for runs whose WGCNA was built before WGCNA_TBmodules_step1.R wrote one.
## The power is not stored in *_WGCNA_workspace.rds, so this script recomputes the scale-free fit table
## from the saved expression matrix (datExpr) with the same rule as WGCNA::pickSoftThreshold
## (unsigned network, powers 1-10 and 12-20, 10 connectivity bins, R^2 cut-off 0.85, fallback power 6
## as in WGCNA_TBmodules_step1.R). It needs only base R.
## Check: for a run that already has the table written by step 1 (TB with HIV), the recomputed table
## must match it (max abs difference < 1e-8), otherwise the script stops.
## Outputs: <run>/10_wgcna_results/<ACC>_WGCNA_soft_threshold.csv (same columns as step 1 writes).
## Each run takes several minutes (one gene-by-gene correlation matrix).

## ---- paths ----
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) return(dirname(normalizePath(ctx$path)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_dir()
ROOT <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, "GSE114192")) && dirname(d) != d) d <- dirname(d); d })

RUNS <- list(c(run = "GSE114192",        acc = "GSE114192"),   # TB vs healthy controls
             c(run = "GSE114192/TB_DM",  acc = "GSE114192"),   # TB with diabetes vs diabetes alone
             c(run = "GSE37250/TB_HIV",  acc = "GSE37250"))    # active vs latent TB with HIV (check run)
POWERS <- c(1:10, seq(12, 20, 2)); R2_CUT <- 0.85; FALLBACK <- 6; N_BREAKS <- 10

## ---- scale-free fit, as in WGCNA::scaleFreeFitIndex ----
scale_free_fit <- function(k, nBreaks = N_BREAKS) {
  discretized.k <- cut(k, nBreaks)
  dk <- tapply(k, discretized.k, mean)
  p.dk <- as.vector(tapply(k, discretized.k, length) / length(k))
  breaks1 <- seq(from = min(k), to = max(k), length = nBreaks + 1)
  hist1 <- hist(k, breaks = breaks1, plot = FALSE, right = TRUE)
  dk2 <- hist1$mids
  dk <- ifelse(is.na(dk), dk2, dk); dk <- ifelse(dk == 0, dk2, dk)
  p.dk <- ifelse(is.na(p.dk), 0, p.dk)
  log.dk <- as.vector(log10(dk)); log.p.dk <- as.numeric(log10(p.dk + 1e-09))
  lm1 <- lm(log.p.dk ~ log.dk); lm2 <- lm(log.p.dk ~ log.dk + I(10^log.dk))
  c(SFT.R.sq = summary(lm1)$r.squared, slope = unname(coef(lm1)[2]), truncated.R.sq = summary(lm2)$adj.r.squared)
}

for (r in RUNS) {
  dir_res <- file.path(ROOT, r[["run"]], "10_wgcna_results")
  ws <- file.path(dir_res, paste0(r[["acc"]], "_WGCNA_workspace.rds"))
  out <- file.path(dir_res, paste0(r[["acc"]], "_WGCNA_soft_threshold.csv"))
  if (!file.exists(ws)) { cat("skip (no workspace):", r[["run"]], "\n"); next }
  datExpr <- as.matrix(readRDS(ws)$datExpr)
  cat(r[["run"]], ":", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")
  a <- abs(cor(datExpr, use = "p"))
  fit <- t(sapply(POWERS, function(b) {
    k <- colSums(a^b) - 1
    c(Power = b, scale_free_fit(k), mean.k. = mean(k), median.k. = median(k), max.k. = max(k))
  }))
  rm(a); gc()
  fit <- as.data.frame(fit)
  ok <- which(fit$SFT.R.sq >= R2_CUT)
  est <- if (length(ok)) POWERS[min(ok)] else NA
  used <- if (is.na(est)) FALLBACK else est
  tab <- data.frame(fit, signed_R2 = -sign(fit$slope) * fit$SFT.R.sq, power_estimate = est,
                    power_used = used, selected = fit$Power == used)
  if (file.exists(out)) {
    old <- read.csv(out)
    num <- c("SFT.R.sq", "slope", "truncated.R.sq", "mean.k.", "median.k.", "max.k.", "signed_R2")
    dev <- max(abs(as.matrix(old[, num]) - as.matrix(tab[, num])))
    same_power <- identical(as.numeric(old$power_used[1]), as.numeric(used))
    cat("  recorded table exists; max abs difference", signif(dev, 3), "; power", used, if (same_power) "(same)" else "(DIFFERENT)", "\n")
    if (dev > 1e-8 || !same_power) stop("Recomputed soft-threshold table does not match the recorded one for ", r[["run"]])
  } else {
    write.csv(tab, out, row.names = FALSE)
    cat("  written:", out, "\n")
  }
  cat("  power estimate:", est, "| power used:", used, "| best fit:", round(max(fit$SFT.R.sq), 3), "\n")
}
