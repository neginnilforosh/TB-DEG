## =============================================================
## WGCNA Pipeline for GSE114192
## =============================================================
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("impute", "preprocessCore", "GO.db", "AnnotationDbi"))
#install.packages("WGCNA")

# Load required libraries
library(WGCNA)
library(dplyr)

# Essential WGCNA settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads() # Enables multi-threading for faster execution

# --- Configuration ---
ACC <- "GSE114192"
BASE_DIR <- file.path("..", ACC)

dir_wgcna   <- file.path(BASE_DIR, "09_wgcna_input")
dir_results <- file.path(BASE_DIR, "10_wgcna_results")
dir_fig     <- file.path(BASE_DIR, "figures")

# Create a new directory for WGCNA numerical outputs
dir.create(dir_results, recursive = TRUE, showWarnings = FALSE)

cat("\n========== STEP 1: LOAD AND FORMAT DATA ==========\n")

# Load Expression Data
expr_file <- file.path(dir_wgcna, paste0(ACC, "_WGCNA_expression_matrix.csv"))
if(!file.exists(expr_file)) stop("Expression matrix not found in 09_wgcna_input folder!")

datExpr <- read.csv(expr_file, row.names = 1, check.names = FALSE)

# WGCNA requires samples as rows and genes as columns.
# If the matrix has genes as rows (typical DESeq2 output), transpose it safely.
if(nrow(datExpr) > ncol(datExpr)) {
    cat("Transposing expression matrix to match WGCNA standards (samples as rows)...\n")
    datExpr <- as.data.frame(t(datExpr))
}
cat("Expression data loaded:", nrow(datExpr), "samples and", ncol(datExpr), "genes.\n")

# Load Trait Data
trait_file <- file.path(dir_wgcna, paste0(ACC, "_WGCNA_trait_file.csv"))
if(!file.exists(trait_file)) stop("Trait file not found in 09_wgcna_input folder!")

datTraits <- read.csv(trait_file, row.names = 1, check.names = FALSE)

# Ensure sample IDs match perfectly between expression matrix and traits
common_samples <- intersect(rownames(datExpr), rownames(datTraits))
if(length(common_samples) == 0) stop("No matching samples found between expression and trait data!")

datExpr <- datExpr[common_samples, ]
datTraits <- datTraits[common_samples, ]
cat("Successfully matched", length(common_samples), "samples between datasets.\n")

cat("\n========== STEP 2: SOFT THRESHOLDING ==========\n")

# Choose a set of candidate soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

cat("Calculating scale-free topology fit indices...\n")
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "unsigned", verbose = 5)

# Plot the Soft Thresholding results
pdf(file.path(dir_fig, paste0(ACC, "_WGCNA_SoftThresholding.pdf")), width = 9, height = 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index plot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n",
     main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red")
abline(h=0.85, col="red") # R-squared threshold line

# Mean connectivity plot
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n",
     main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

# Automatically select the optimal power
softPower <- sft$powerEstimate
if(is.na(softPower)) {
    cat("WARNING: Algorithm could not determine an optimal soft power. Defaulting to 6.\n")
    softPower <- 6
}
cat(">>> Selected Soft Threshold Power:", softPower, "\n")

cat("\n========== STEP 3: NETWORK CONSTRUCTION ==========\n")
cat("Constructing gene network and identifying modules (this may take a few minutes)...\n")

# One-step network construction and module detection
net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE,
                        verbose = 3)

# Convert numeric module labels into WGCNA standard colors
moduleColors <- labels2colors(net$colors)
cat(">>> Found", length(table(moduleColors)), "modules.\n")

# Plot the dendrogram with module colors
pdf(file.path(dir_fig, paste0(ACC, "_WGCNA_ModuleDendrogram.pdf")), width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste("Gene Dendrogram and Module Colors (", ACC, ")", sep=""))
dev.off()

cat("\n========== STEP 4: MODULE-TRAIT CORRELATION ==========\n")

nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# Recalculate Module Eigengenes (MEs) with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

cat("Calculating Pearson correlations between modules and clinical traits...\n")
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Create text matrix for the heatmap
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Plot Module-Trait Correlation Heatmap
pdf(file.path(dir_fig, paste0(ACC, "_WGCNA_ModuleTrait_Heatmap.pdf")), width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships (", ACC, ")", sep=""))
dev.off()

cat("\n========== STEP 5: EXPORT RESULTS ==========\n")

# Export module assignments for all genes
gene_module_info <- data.frame(
    Gene = colnames(datExpr),
    ModuleColor = moduleColors
)
write.csv(gene_module_info, file.path(dir_results, paste0(ACC, "_Gene_Module_Assignment.csv")), row.names = FALSE)

