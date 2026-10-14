## Reproducibility record: software/package versions, important parameter settings, thresholds,
## database versions where available, and random seeds.
##
## Reads every script under scripts/ and writes, into reproducibility/ at the repository root:
##   package_versions.csv   every package the scripts use, with the installed version (or "not installed")
##   sessionInfo.txt        R version, platform, BLAS/LAPACK, locale
##   PARAMETERS.md          every script's CONFIG block, the run presets in config.R, all set.seed() calls,
##                          and annotation / database versions
## Re-run it whenever packages are updated or a CONFIG value changes, and commit the folder.
## Run from anywhere inside the repository:   Rscript scripts/record_reproducibility.R

here <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
ROOT <- local({ d <- here; while (!dir.exists(file.path(d, "scripts")) && dirname(d) != d) d <- dirname(d); d })
scripts <- sort(list.files(file.path(ROOT, "scripts"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE))
scripts <- scripts[!grepl("record_reproducibility|remerge_test", scripts)]
out <- file.path(ROOT, "reproducibility"); dir.create(out, showWarnings = FALSE)
rel <- function(p) sub(paste0("^", ROOT, "/"), "", p)
cat(">>>", length(scripts), "scripts scanned\n")

## ---- 1. packages ----
src <- unlist(lapply(scripts, readLines, warn = FALSE))
src <- sub("#.*$", "", src)                                           # ignore comments
pk <- c(unlist(regmatches(src, gregexpr("(?<=\\blibrary\\()[A-Za-z0-9.]+", src, perl = TRUE))),
        unlist(regmatches(src, gregexpr("(?<=\\brequire\\()[A-Za-z0-9.]+", src, perl = TRUE))),
        unlist(regmatches(src, gregexpr("[A-Za-z][A-Za-z0-9.]+(?=::)", src, perl = TRUE))),
        unlist(regmatches(src, gregexpr("(?<=requireNamespace\\(\")[A-Za-z0-9.]+", src, perl = TRUE))))
## packages named inside c("...") vectors of required packages
for (m in regmatches(src, gregexpr("(req|required_pkgs|pkgs?)\\s*<-\\s*c\\([^)]*\\)", src)))
  for (x in m) pk <- c(pk, regmatches(x, gregexpr("(?<=\")[A-Za-z][A-Za-z0-9.]+(?=\")", x, perl = TRUE))[[1]])
base_pk <- rownames(installed.packages(priority = "base"))
pk <- sort(setdiff(unique(pk), c(base_pk, "utils", "stats", "grDevices", "graphics", "methods", "tools")))
ver <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"), "")
pv <- data.frame(Package = pk, Version = unname(ver), stringsAsFactors = FALSE)
write.csv(pv, file.path(out, "package_versions.csv"), row.names = FALSE)
cat(">>>", nrow(pv), "packages recorded;", sum(pv$Version == "not installed"), "not installed here\n")

## ---- 2. session ----
writeLines(capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))

## ---- 3. parameters ----
md <- c("# Reproducibility record", "", paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " with ",
        R.version.string, " on ", utils::sessionInfo()$running, "."), "",
        "Package versions: `package_versions.csv`. Full session: `sessionInfo.txt`.", "")

## annotation / database versions that can be read from installed packages
db <- c(org.Hs.eg.db = "org.Hs.eg.db", msigdbr = "msigdbr", msigdbdf = "msigdbdf", STRINGdb = "STRINGdb",
        clusterProfiler = "clusterProfiler", ChIPseeker = "ChIPseeker",
        TxDb_hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene", drugfindR = "drugfindR", WGCNA = "WGCNA", DESeq2 = "DESeq2")
md <- c(md, "## Annotation and database versions", "", "| Resource | Version |", "|---|---|")
for (k in names(db)) md <- c(md, sprintf("| %s | %s |", db[[k]],
  tryCatch(as.character(utils::packageVersion(db[[k]])), error = function(e) "not installed")))
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  meta <- tryCatch(AnnotationDbi::metadata(org.Hs.eg.db::org.Hs.eg.db), error = function(e) NULL)
  if (!is.null(meta)) { keep <- meta[meta$name %in% c("EGSOURCEDATE", "ENSOURCEDATE", "GOSOURCEDATE"), ]
    md <- c(md, "", paste0("org.Hs.eg.db source dates: ", paste(keep$name, keep$value, sep = " = ", collapse = "; "))) }
}
md <- c(md, "", "Web resources queried at run time (versions are those live on the run date): STRING, iLINCS,",
        "Enrichr, DGIdb, KEGG REST, Metabolomics Workbench. See the run dates in each output folder.", "")

## presets
cfg <- file.path(ROOT, "scripts", "wgcna_network", "config.R")
if (file.exists(cfg)) {
  x <- readLines(cfg, warn = FALSE); a <- grep("^PRESETS", x); b <- a
  while (b < length(x) && !grepl("^\\)", x[b])) b <- b + 1
  md <- c(md, "## Run presets (scripts/wgcna_network/config.R)", "", "```r", x[a:b], "```", "")
}

## seeds
seeds <- do.call(rbind, lapply(scripts, function(f) { x <- readLines(f, warn = FALSE); i <- grep("set\\.seed\\(", x)
  if (length(i)) data.frame(Script = rel(f), Line = i, Call = trimws(x[i])) else NULL }))
md <- c(md, "## Random seeds", "", "| Script | Line | Call |", "|---|---|---|",
        if (!is.null(seeds)) sprintf("| %s | %d | `%s` |", seeds$Script, seeds$Line, seeds$Call) else "| none | | |", "")

## every CONFIG block
md <- c(md, "## Parameters by script", "")
for (f in scripts) {
  x <- readLines(f, warn = FALSE); a <- grep("^## ---- CONFIG", x)
  if (!length(a)) next
  b <- a + 1; while (b <= length(x) && !grepl("^## ----", x[b])) b <- b + 1
  blk <- x[(a + 1):(b - 1)]; blk <- blk[nzchar(trimws(blk))]
  if (length(blk)) md <- c(md, paste0("### ", rel(f)), "", "```r", blk, "```", "")
}
writeLines(md, file.path(out, "PARAMETERS.md"))
cat(">>> wrote", rel(file.path(out, "PARAMETERS.md")), "| seeds found:", if (is.null(seeds)) 0 else nrow(seeds), "\n")
cat(">>> DONE. Commit the reproducibility/ folder.\n")
