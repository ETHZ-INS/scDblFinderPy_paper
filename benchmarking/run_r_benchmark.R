suppressMessages({
  library(SingleCellExperiment)
  library(scds)
  library(scDblFinder)
  library(Seurat)
  library(DoubletFinder)
})

ALL_METHODS <- c("scDblFinder.clusters", "scDblFinder.random",
                  "bcds", "cxds", "hybrid", "computeDoubletDensity",
                  "DoubletFinder")

# Which methods to actually (re)compute; everything else is carried over
# from this script's last run (see below). Usage:
#   Rscript run_r_benchmark.R                     # DoubletFinder only (default)
#   Rscript run_r_benchmark.R --all                # every method
#   Rscript run_r_benchmark.R bcds,cxds,hybrid      # comma-separated subset
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
  methods_to_run <- c("DoubletFinder")
} else if (identical(args[1], "--all")) {
  methods_to_run <- ALL_METHODS
} else {
  methods_to_run <- strsplit(args[1], ",")[[1]]
}
unknown <- setdiff(methods_to_run, ALL_METHODS)
if (length(unknown) > 0) {
  stop("Unknown method(s): ", paste(unknown, collapse=", "),
       "\nValid methods: ", paste(ALL_METHODS, collapse=", "))
}
cat("Methods to (re)run:", paste(methods_to_run, collapse=", "), "\n")

# benchmark.results_R.rds is self-perpetuating: each run reads whatever it
# last wrote (never any static original-paper file, which goes stale the
# moment any method's package version moves on) and carries forward the
# methods it isn't asked to recompute this time. The very first run should
# pass --all so every method has a real, current value to carry forward;
# otherwise the methods you skip simply won't be in the output yet.
result_file <- "benchmark.results_R.rds"
if (file.exists(result_file)) {
  cat("Loading previous R-package benchmark results from", result_file, "...\n")
  old <- readRDS(result_file)
} else {
  if (!setequal(methods_to_run, ALL_METHODS)) {
    cat("WARNING:", result_file, "doesn't exist yet and you're not rerunning",
        "every method (--all). The output will be missing rows for:",
        paste(setdiff(ALL_METHODS, methods_to_run), collapse=", "), "\n")
  }
  old <- data.frame(dataset=character(), method=character(),
                     AUPRC=numeric(), AUROC=numeric(), elapsed=numeric())
}

# Chord is always dropped (unmaintained); Scrublet and Vaeda moved to
# run_scrublet_benchmark.py/run_vaeda_benchmark.py (they're Python tools).
unchanged_methods <- setdiff(ALL_METHODS, methods_to_run)
kept <- old[old$method %in% unchanged_methods, ]

rds_files <- list.files(file.path("datasets"), pattern="\\.rds$", full.names=TRUE)

# scDblFinder(): clusters=TRUE auto-clusters and generates artificial
# doublets between clusters; clusters=NULL generates purely random
# artificial doublets. This is exactly the clusters/random distinction
# used throughout this benchmark (matching scDblFinder.Py.clusters/random).
run_scdblfinder <- function(counts, clusters) {
  sce <- SingleCellExperiment(assays=list(counts=counts))
  res <- scDblFinder(sce, clusters=clusters, verbose=FALSE)
  colData(res)$scDblFinder.score
}

method_fns <- list(
  "scDblFinder.clusters" = function(counts) run_scdblfinder(counts, clusters=TRUE),
  "scDblFinder.random"   = function(counts) run_scdblfinder(counts, clusters=NULL),

  "bcds" = function(counts) {
    sce <- SingleCellExperiment(assays=list(counts=counts))
    sce <- scds::bcds(sce, verb=FALSE)
    colData(sce)$bcds_score
  },
  "cxds" = function(counts) {
    sce <- SingleCellExperiment(assays=list(counts=counts))
    sce <- scds::cxds(sce, verb=FALSE)
    colData(sce)$cxds_score
  },
  "hybrid" = function(counts) {
    sce <- SingleCellExperiment(assays=list(counts=counts))
    sce <- scds::cxds_bcds_hybrid(sce, verb=FALSE)
    colData(sce)$hybrid_score
  },
  "computeDoubletDensity" = function(counts) {
    scDblFinder::computeDoubletDensity(counts)
  },

  "DoubletFinder" = function(counts) {
    seu <- CreateSeuratObject(counts=counts)
    seu <- NormalizeData(seu, verbose=FALSE)
    seu <- FindVariableFeatures(seu, verbose=FALSE)
    seu <- ScaleData(seu, verbose=FALSE)
    seu <- RunPCA(seu, npcs=30, verbose=FALSE)
    seu <- FindNeighbors(seu, dims=1:10, verbose=FALSE)
    seu <- FindClusters(seu, resolution=0.8, verbose=FALSE)

    sweep.res <- paramSweep(seu, PCs=1:10, sct=FALSE)
    sweep.stats <- summarizeSweep(sweep.res, GT=FALSE)
    # find.pK() unconditionally calls plot() with no way to suppress it;
    # under Rscript that falls back to the default pdf device and dumps
    # an Rplots.pdf into the cwd unless we route it to a null device.
    pdf(file=nullfile())
    pk.stats <- find.pK(sweep.stats)
    dev.off()
    best.pK <- as.numeric(as.character(pk.stats$pK[which.max(pk.stats$BCmetric)]))

    # 10x-style heuristic (~0.8% doublets per 1000 cells) rather than the
    # known ground truth, since real experiments don't have that available.
    homotypic.prop <- modelHomotypic(seu@meta.data$seurat_clusters)
    doublet_rate <- 0.008 * (ncol(seu) / 1000)
    nExp_poi <- round(doublet_rate * ncol(seu))
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

    # reuse.pANN must be omitted, not set to FALSE: DoubletFinder 2.0.6
    # checks `!is.null(reuse.pANN)`, and is.null(FALSE) is FALSE, so passing
    # FALSE wrongly takes the "reuse old pANN" branch and errors out.
    seu <- doubletFinder(seu, PCs=1:10, pN=0.25, pK=best.pK, nExp=nExp_poi.adj, sct=FALSE)
    pann_col <- grep("^pANN", colnames(seu@meta.data), value=TRUE)
    seu@meta.data[[pann_col]]
  }
)

run_method <- function(method, counts, ds_name) {
  st <- Sys.time()
  scores <- tryCatch(method_fns[[method]](counts), error = function(e) {
    cat("    Failed on", method, "for", ds_name, ":", conditionMessage(e), "\n")
    NULL
  })
  list(scores=scores, elapsed=as.numeric(Sys.time() - st, units="secs"))
}

new_results <- list()
for (f in rds_files) {
  ds_name <- sub("\\.rds$", "", basename(f))

  obj <- readRDS(f)
  counts <- obj[[1]]
  truth_labels <- as.integer(tolower(obj[[2]]) == "doublet")

  for (method in methods_to_run) {
    cat("Evaluating", ds_name, "-", method, "...\n")
    res <- run_method(method, counts, ds_name)
    if (!is.null(res$scores)) {
      auroc <- as.numeric(pROC::auc(pROC::roc(truth_labels, res$scores, quiet=TRUE)))
      pr <- PRROC::pr.curve(scores.class0=res$scores, weights.class0=truth_labels)
      auprc <- pr$auc.integral

      new_results[[length(new_results) + 1]] <- data.frame(
        dataset=ds_name, method=method,
        AUPRC=auprc, AUROC=auroc, elapsed=res$elapsed
      )
    }
  }
}

new_rows <- do.call(rbind, new_results)
e <- rbind(kept, new_rows)

saveRDS(e, result_file)
cat("Saved R-package benchmark results to", result_file, "\n")
