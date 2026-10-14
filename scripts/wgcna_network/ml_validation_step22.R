## STEP 22: machine-learning validation (Ratul, 7 Sep: "you can add the ML validation step after this").
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

## ---- CONFIG ----
TIERS     <- c("High")              # genes used as features; c("High","Medium") for a larger panel
TESTS <- list(                       # dataset, group A (control), group B (case), label
  list(ds = "GSE161829", a = "TBneg",   b = "ATB",    lab = "GSE161829 ATB vs TBneg"),
  list(ds = "GSE161829", a = "LTBI",    b = "ATB",    lab = "GSE161829 ATB vs LTBI"),
  list(ds = "GSE161829", a = "TBneg",   b = "LTBI",   lab = "GSE161829 LTBI vs TBneg"),
  list(ds = "GSE222001", a = "Healthy", b = "Active", lab = "GSE222001 Active vs Healthy (n=5+5)"),
  list(ds = "GSE222001", a = "Latent",  b = "Active", lab = "GSE222001 Active vs Latent (n=5+5)"))
SEED <- 42; N_BOOT <- 1000; N_RANDOM <- 1000; ALPHA <- 0.5

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
add(paste0("Discovery (", ACC, ") ", CASE, " vs ", CONTROL, " [in-sample: genes chosen here]"), "Signature score", sig_score(Dz, up, dn), Dy,
    sapply(rand_sets, function(g) { z <- zscore(Dx[g, ]); s <- colMeans(z); max(auc(s, Dy), auc(-s, Dy)) }))
add(paste0("Discovery (", ACC, ") ", CASE, " vs ", CONTROL, " [in-sample: genes chosen here]"), "Elastic net (5-fold CV)", cv_pred, Dy)

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
