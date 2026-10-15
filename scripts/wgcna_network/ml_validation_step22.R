## STEP 22: machine-learning validation.
##
## Question: do the prioritised genes, learned in the discovery cohort, separate TB from controls in
## cohorts they have never seen -- and do they separate ACTIVE from LATENT TB?
## Two classifiers on the same genes:
##   A. untrained signature score = mean z(up genes) - mean z(down genes); no fitting at all, so it
##      cannot overfit -- the direction of each gene comes from the discovery DEG table
##   B. elastic-net logistic regression (glmnet, alpha 0.5) trained on the discovery cohort only
## Each external dataset is z-scored gene by gene on its own samples, so platform/scale differences
## between studies do not leak into the comparison. AUCs carry bootstrap 95% CIs, and the score's AUC is
## compared with 1,000 random gene sets of the same size, so "AUC 0.9" is read against what any blood
## gene set achieves. Outputs in <run>/24_ml_validation/.
##
## Leakage (D12): the prioritised genes were chosen on ALL discovery samples, so their discovery AUC is
## in-sample and optimistic; it is kept only as a descriptive row. The honest discovery estimate is a
## repeated nested cross-validation of the whole signature-building procedure: in every outer training
## fold the genes are re-selected (Welch t-test, same panel size), standardised with training-fold
## parameters, the elastic-net lambda is chosen by inner CV, and only then is the held-out fold predicted.
## External datasets are true hold-outs: nothing is fitted on them.

## ---- CONFIG ----
TIERS     <- c("High")              # genes used as features; c("High","Medium") for a larger panel
TESTS <- list(                       # dataset, group A (control), group B (case), label
  list(ds = "GSE161829", a = "TBneg",   b = "ATB",    lab = "GSE161829 ATB vs TBneg"),
  list(ds = "GSE161829", a = "LTBI",    b = "ATB",    lab = "GSE161829 ATB vs LTBI"),
  list(ds = "GSE161829", a = "TBneg",   b = "LTBI",   lab = "GSE161829 LTBI vs TBneg"),
  list(ds = "GSE222001", a = "Healthy", b = "Active", lab = "GSE222001 Active vs Healthy (n=5+5)"),
  list(ds = "GSE222001", a = "Latent",  b = "Active", lab = "GSE222001 Active vs Latent (n=5+5)"))
SEED <- 42; N_BOOT <- 1000; N_RANDOM <- 1000; ALPHA <- 0.5
N_OUTER <- 5; N_REPEATS <- 10        # nested CV: 5 outer folds, repeated 10 times

for (p in c("glmnet", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
suppressPackageStartupMessages(library(ggplot2))
get_script_dir <- function() { a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() }
SCRIPT_DIR <- get_script_dir()
source(local({ d <- SCRIPT_DIR; while (!file.exists(file.path(d, "config.R")) && dirname(d) != d) d <- dirname(d); file.path(d, "config.R") }))
ROOT_DIR <- local({ d <- SCRIPT_DIR; while (!dir.exists(file.path(d, ACC, "02_metadata")) && dirname(d) != d) d <- dirname(d); d })
BASE_DIR <- file.path(ROOT_DIR, ACC); OUT_DIR <- run_out_dir(BASE_DIR)
dir_out  <- file.path(OUT_DIR, "24_ml_validation"); dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
FIG_DPI <- 600
save_gg <- function(f, p, w, h) { ggsave(f, p, width = w, height = h, dpi = FIG_DPI); ggsave(sub("\\.png$", ".pdf", f), p, width = w, height = h) }

## ---- data helpers ----
ens_key <- function(id) sub("\\..*$", "", sub("\\|.*$", "", id))
load_dataset <- function(base) {
  vf <- list.files(file.path(base, "04_normalized"), pattern = "_vst_normalized_matrix\\.csv$", full.names = TRUE)[1]
  mf <- list.files(file.path(base, "02_metadata"), pattern = "_sample_metadata\\.csv$", full.names = TRUE)[1]
  if (is.na(vf) || is.na(mf)) stop("No normalized matrix / metadata in ", base)
  x <- as.matrix(read.csv(vf, row.names = 1, check.names = FALSE)); md <- read.csv(mf, stringsAsFactors = FALSE)
  k <- ens_key(rownames(x)); keep <- grepl("^ENSG", k) & !duplicated(k); x <- x[keep, , drop = FALSE]; rownames(x) <- k[keep]
  grp <- setNames(as.character(md$group), md$sample_id); x <- x[, intersect(colnames(x), names(grp)), drop = FALSE]
  list(x = x, grp = grp[colnames(x)])
}
zscore <- function(x) { z <- t(scale(t(x))); z[is.na(z)] <- 0; z }
auc <- function(score, y) { r <- rank(score); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
boot_ci <- function(score, y) { set.seed(SEED); i1 <- which(y == 1); i0 <- which(y == 0)
  b <- replicate(N_BOOT, { s <- c(sample(i1, replace = TRUE), sample(i0, replace = TRUE)); auc(score[s], y[s]) })
  stats::quantile(b, c(0.025, 0.975), na.rm = TRUE) }
sig_score <- function(z, up, dn) { s <- 0
  if (length(up)) s <- s + colMeans(z[up, , drop = FALSE]); if (length(dn)) s <- s - colMeans(z[dn, , drop = FALSE]); s }

## ---- discovery ----
pri <- read.csv(file.path(OUT_DIR, "19_consensus", paste0(ACC, "_GenePriority_all_modules.csv")), stringsAsFactors = FALSE)
feat <- pri[pri$Final_Priority %in% TIERS & !is.na(pri$log2FC), ]
D <- load_dataset(in_dir(BASE_DIR, "02_metadata") |> dirname())
D_keep <- D$grp %in% c(CONTROL, CASE); Dx <- D$x[, D_keep]; Dy <- as.integer(D$grp[D_keep] == CASE)
ext <- lapply(unique(sapply(TESTS, `[[`, "ds")), function(ds) load_dataset(file.path(ROOT_DIR, ds)))
names(ext) <- unique(sapply(TESTS, `[[`, "ds"))
common <- Reduce(intersect, c(list(rownames(Dx)), lapply(ext, function(e) rownames(e$x))))
genes <- intersect(feat$Ensembl, common)
cat(">>> run", run_id, "|", CONTROL, "(", sum(Dy == 0), ") vs", CASE, "(", sum(Dy == 1), ") |",
    nrow(feat), paste(TIERS, collapse = "+"), "genes,", length(genes), "measured in every test dataset\n")
if (length(genes) < 5) stop("Fewer than 5 feature genes are measured in all datasets.")
up <- genes[genes %in% feat$Ensembl[feat$log2FC > 0]]; dn <- setdiff(genes, up)
Dz <- zscore(Dx[genes, , drop = FALSE])

## elastic net on the discovery cohort only; CV predictions give an honest discovery AUC
set.seed(SEED)
cv <- glmnet::cv.glmnet(t(Dz), Dy, family = "binomial", alpha = ALPHA, nfolds = 5, keep = TRUE, type.measure = "auc")
lam_i <- which(cv$lambda == cv$lambda.min)
cv_pred <- cv$fit.preval[, lam_i]
coefs <- as.matrix(stats::coef(cv, s = "lambda.min")); coefs <- coefs[coefs[, 1] != 0, , drop = FALSE]
cf <- data.frame(Ensembl = rownames(coefs), Coefficient = round(coefs[, 1], 4))
cf$Gene <- pri$Gene[match(cf$Ensembl, pri$Ensembl)]
write.csv(cf, file.path(dir_out, paste0(ACC, "_elasticnet_coefficients.csv")), row.names = FALSE)
cat(">>> elastic net kept", nrow(cf) - 1, "of", length(genes), "genes\n")

## random gene sets for the score's null: same size, drawn from all genes measured everywhere
set.seed(SEED); pool <- setdiff(common, genes)
rand_sets <- replicate(N_RANDOM, sample(pool, length(genes)), simplify = FALSE)

res <- list(); roc <- list()
add <- function(label, model, score, y, null = NULL) {
  a <- auc(score, y); ci <- boot_ci(score, y)
  res[[length(res) + 1]] <<- data.frame(Test = label, Model = model, N_case = sum(y == 1), N_control = sum(y == 0),
    AUC = round(a, 3), CI_low = round(ci[1], 3), CI_high = round(ci[2], 3),
    Random_gene_set_median_AUC = if (is.null(null)) NA else round(stats::median(null), 3),
    P_vs_random = if (is.null(null)) NA else signif((sum(null >= a) + 1) / (length(null) + 1), 3))
  o <- order(-score); tp <- cumsum(y[o] == 1) / sum(y == 1); fp <- cumsum(y[o] == 0) / sum(y == 0)
  roc[[length(roc) + 1]] <<- data.frame(Test = label, Model = model, FPR = c(0, fp), TPR = c(0, tp))
}
## discovery
add(paste0("Discovery (", ACC, ") ", CASE, " vs ", CONTROL, " [in-sample, optimistic: genes chosen on these samples]"), "Signature score", sig_score(Dz, up, dn), Dy,
    sapply(rand_sets, function(g) { z <- zscore(Dx[g, ]); s <- colMeans(z); max(auc(s, Dy), auc(-s, Dy)) }))
add(paste0("Discovery (", ACC, ") ", CASE, " vs ", CONTROL, " [in-sample, optimistic: genes chosen on these samples]"), "Elastic net (5-fold CV)", cv_pred, Dy)

## ---- nested cross-validation of the whole signature-building procedure (leakage-free discovery AUC) ----
welch <- function(x, y) {                     # row-wise Welch t statistic and two-sided p-value
  a <- x[, y == 1, drop = FALSE]; b <- x[, y == 0, drop = FALSE]
  na <- ncol(a); nb <- ncol(b); va <- apply(a, 1, stats::var); vb <- apply(b, 1, stats::var)
  se <- sqrt(va / na + vb / nb); tt <- (rowMeans(a) - rowMeans(b)) / se
  df <- se^4 / ((va / na)^2 / (na - 1) + (vb / nb)^2 / (nb - 1))
  pp <- 2 * stats::pt(-abs(tt), df)
  tt[!is.finite(tt)] <- 0; pp[!is.finite(pp)] <- 1
  names(tt) <- rownames(x); names(pp) <- rownames(x)
  list(t = tt, p = pp)
}
strat_folds <- function(y, k) {               # stratified fold labels 1..k
  f <- integer(length(y))
  for (cl in c(0, 1)) { i <- which(y == cl); f[i] <- sample(rep_len(seq_len(k), length(i))) }
  f
}
k_feat <- length(genes)                       # same panel size as the prioritised genes
Xn <- Dx[common, , drop = FALSE]              # candidate genes: everything measured in every dataset
set.seed(SEED)
nest <- list()
for (r in seq_len(N_REPEATS)) {
  folds <- strat_folds(Dy, N_OUTER); p_sig <- p_en <- rep(NA_real_, length(Dy))
  for (k in seq_len(N_OUTER)) {
    tr <- folds != k; te <- folds == k
    w <- welch(Xn[, tr, drop = FALSE], Dy[tr])
    sel <- names(sort(w$p))[seq_len(min(k_feat, length(w$p)))]
    mu <- rowMeans(Xn[sel, tr, drop = FALSE]); sdv <- apply(Xn[sel, tr, drop = FALSE], 1, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    ztr <- (Xn[sel, tr, drop = FALSE] - mu) / sdv
    zte <- (Xn[sel, te, drop = FALSE] - mu) / sdv
    up_k <- sel[w$t[sel] > 0]; dn_k <- setdiff(sel, up_k)
    p_sig[te] <- sig_score(zte, up_k, dn_k)
    fit <- glmnet::cv.glmnet(t(ztr), Dy[tr], family = "binomial", alpha = ALPHA, nfolds = 5, type.measure = "auc")
    p_en[te] <- as.numeric(stats::predict(fit, newx = t(zte), s = "lambda.min", type = "link"))
  }
  nest[[r]] <- data.frame(Repeat = r, AUC_signature = auc(p_sig, Dy), AUC_elastic_net = auc(p_en, Dy))
}
NEST <- do.call(rbind, nest)
write.csv(NEST, file.path(dir_out, paste0(ACC, "_ml_nested_cv.csv")), row.names = FALSE)
nest_lab <- paste0("Discovery (", ACC, ") ", CASE, " vs ", CONTROL, " [nested CV ", N_OUTER, "-fold x ", N_REPEATS,
                   "; genes re-selected in each training fold; CI = 2.5-97.5% of repeats]")
for (m in c("AUC_signature", "AUC_elastic_net")) {
  v <- NEST[[m]]
  res[[length(res) + 1]] <- data.frame(Test = nest_lab, Model = if (m == "AUC_signature") "Signature score (nested CV)" else "Elastic net (nested CV)",
    N_case = sum(Dy == 1), N_control = sum(Dy == 0), AUC = round(mean(v), 3),
    CI_low = round(unname(stats::quantile(v, 0.025)), 3), CI_high = round(unname(stats::quantile(v, 0.975)), 3),
    Random_gene_set_median_AUC = NA, P_vs_random = NA)
}
cat(">>> nested CV (", N_OUTER, "x", N_REPEATS, "): signature AUC", round(mean(NEST$AUC_signature), 3),
    "| elastic net AUC", round(mean(NEST$AUC_elastic_net), 3), "\n")

## external tests
for (t in TESTS) {
  e <- ext[[t$ds]]; keep <- e$grp %in% c(t$a, t$b)
  if (sum(e$grp == t$a) < 3 || sum(e$grp == t$b) < 3) { cat("!! too few samples:", t$lab, "\n"); next }
  x <- e$x[, keep]; y <- as.integer(e$grp[keep] == t$b); z <- zscore(x[genes, , drop = FALSE])
  ## random sets: the random set's direction is unknown, so take its better orientation (a generous null)
  null <- sapply(rand_sets, function(g) { s <- colMeans(zscore(x[g, , drop = FALSE])); max(auc(s, y), auc(-s, y)) })
  add(t$lab, "Signature score", sig_score(z, up, dn), y, null)
  add(t$lab, "Elastic net (trained on discovery)", as.numeric(stats::predict(cv, newx = t(z), s = "lambda.min", type = "link")), y)
}
R <- do.call(rbind, res); write.csv(R, file.path(dir_out, paste0(ACC, "_ml_validation_AUC.csv")), row.names = FALSE)
cat("\n========== AUC ==========\n"); print(R, row.names = FALSE)

RC <- do.call(rbind, roc); RC$Test <- factor(RC$Test, levels = unique(RC$Test))
p <- ggplot(RC, aes(FPR, TPR, colour = Model)) + geom_abline(linetype = 2, colour = "grey60") + geom_step() +
  facet_wrap(~ Test, ncol = 3) + coord_equal() + theme_bw(base_size = 8) + theme(legend.position = "bottom") +
  labs(title = paste0(ACC, ": ", length(genes), " prioritised genes, trained on discovery only"), x = "False positive rate", y = "True positive rate")
save_gg(file.path(dir_out, paste0(ACC, "_ml_validation_ROC.png")), p, 8, 6.5)
cat("\nNote: GSE222001 has 5 samples per group, so its AUCs and CIs are rough. Random gene sets are scored in their\n",
    "better orientation, which favours the null -- an AUC above them is conservative evidence.\n", sep = "")
cat(">>> DONE. Outputs in", dir_out, "\n")
