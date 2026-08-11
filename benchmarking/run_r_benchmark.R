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

# Args:
#   Rscript run_r_benchmark.R                                  # DoubletFinder only, seeds 42,0,1 (default)
#   Rscript run_r_benchmark.R --all                             # every method, seeds 42,0,1
#   Rscript run_r_benchmark.R bcds,cxds,hybrid                  # comma-separated subset
#   Rscript run_r_benchmark.R --all --seeds 42,0,1              # explicit seeds
# --seeds (if present) is stripped out before the methods-selection arg is parsed.
args <- commandArgs(trailingOnly=TRUE)
seeds_idx <- which(args == "--seeds")
if (length(seeds_idx) > 0) {
  SEEDS <- as.integer(strsplit(args[seeds_idx + 1], ",")[[1]])
  args <- args[-c(seeds_idx, seeds_idx + 1)]
} else {
  SEEDS <- c(42L, 0L, 1L)
}

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
cat("Seeds:", paste(SEEDS, collapse=", "), "\n")

# benchmark.results_R.rds is self-perpetuating: each run reads whatever it
# last wrote (never any static original-paper file, which goes stale the
# moment any method's package version moves on) and carries forward the
# methods it isn't asked to recompute this time. The very first run should
# pass --all so every method has a real, current value to carry forward;
# otherwise the methods you skip simply won't be in the output yet.
# benchmark.percell_scores_R.rds mirrors the same carry-over scheme, one row
# per (dataset, method, seed, cell) - this is what lets a later analysis ask
# whether cell-wise averaging the score across the 3 seeds beats the average
# of the 3 seeds' individual AUPRCs (Pierre-Luc's question).
result_file <- "benchmark.results_R.rds"
percell_file <- "benchmark.percell_scores_R.rds"

empty_results <- data.frame(dataset=character(), method=character(), seed=integer(),
                             AUPRC=numeric(), AUROC=numeric(), elapsed=numeric())
empty_percell <- data.frame(dataset=character(), method=character(), seed=integer(),
                             cell=character(), score=numeric(), truth=integer())

load_cache <- function(path, empty, label) {
  if (!file.exists(path)) {
    if (!setequal(methods_to_run, ALL_METHODS)) {
      cat("WARNING:", path, "doesn't exist yet and you're not rerunning",
          "every method (--all). The output will be missing rows for:",
          paste(setdiff(ALL_METHODS, methods_to_run), collapse=", "), "\n")
    }
    return(empty)
  }
  cat("Loading previous", label, "from", path, "...\n")
  old <- readRDS(path)
  if (!"seed" %in% colnames(old)) {
    cat("WARNING:", path, "predates per-seed tracking (no 'seed' column) -",
        "discarding it. Rerun with --all to repopulate every method.\n")
    return(empty)
  }
  old
}

old_results <- load_cache(result_file, empty_results, "R-package benchmark results")
old_percell <- load_cache(percell_file, empty_percell, "per-cell doublet scores")

# Chord is always dropped (unmaintained); Scrublet and Vaeda moved to
# run_scrublet_benchmark.py/run_vaeda_benchmark.py (they're Python tools).
unchanged_methods <- setdiff(ALL_METHODS, methods_to_run)
kept_results <- old_results[old_results$method %in% unchanged_methods, ]
kept_percell <- old_percell[old_percell$method %in% unchanged_methods, ]

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

# method_fns take (counts, seed) - every method sets a local seed itself
# (rather than relying solely on a global set.seed() before the call) since
# a couple of Seurat steps inside DoubletFinder default seed.use/random.seed
# to fixed constants (42 / 0) that ignore the global RNG state entirely.
method_fns <- list(
  "scDblFinder.clusters" = function(counts, seed) { set.seed(seed); run_scdblfinder(counts, clusters=TRUE) },
  "scDblFinder.random"   = function(counts, seed) { set.seed(seed); run_scdblfinder(counts, clusters=NULL) },

  "bcds" = function(counts, seed) {
    set.seed(seed)
    sce <- SingleCellExperiment(assays=list(counts=counts))
    sce <- scds::bcds(sce, verb=FALSE)
    colData(sce)$bcds_score
  },
  "cxds" = function(counts, seed) {
    set.seed(seed)
    sce <- SingleCellExperiment(assays=list(counts=counts))
    sce <- scds::cxds(sce, verb=FALSE)
    colData(sce)$cxds_score
  },
  "hybrid" = function(counts, seed) {
    set.seed(seed)
    sce <- SingleCellExperiment(assays=list(counts=counts))
    sce <- scds::cxds_bcds_hybrid(sce, verb=FALSE)
    colData(sce)$hybrid_score
  },
  "computeDoubletDensity" = function(counts, seed) {
    set.seed(seed)
    scDblFinder::computeDoubletDensity(counts)
  },

  "DoubletFinder" = function(counts, seed) {
    set.seed(seed)
    seu <- CreateSeuratObject(counts=counts)
    seu <- NormalizeData(seu, verbose=FALSE)
    seu <- FindVariableFeatures(seu, verbose=FALSE)
    seu <- ScaleData(seu, verbose=FALSE)
    seu <- RunPCA(seu, npcs=30, verbose=FALSE, seed.use=seed)
    seu <- FindNeighbors(seu, dims=1:10, verbose=FALSE)
    seu <- FindClusters(seu, resolution=0.8, verbose=FALSE, random.seed=seed)

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

    # paramSweep()'s internal parallel_paramSweep() rebuilds its own Seurat
    # object per (pN, pK) combo and calls RunPCA() on it WITHOUT seed.use,
    # so it always resets R's global RNG via Seurat's hardcoded seed.use=42
    # default - clobbering whatever `seed` set beforehand, regardless of its
    # value. Nothing between paramSweep() finishing and here touches the
    # RNG again (summarizeSweep/find.pK/modelHomotypic don't seed), so
    # without reseeding here, doubletFinder()'s own sample() calls for
    # artificial-doublet generation would be identical across all 3 seeds.
    # (best.pK selection above is still effectively seed-independent either
    # way, since it's entirely downstream of that clobbered paramSweep RNG
    # state - only the final doublet scores computed below are recoverably
    # seeded.)
    set.seed(seed)
    # reuse.pANN must be omitted, not set to FALSE: DoubletFinder 2.0.6
    # checks `!is.null(reuse.pANN)`, and is.null(FALSE) is FALSE, so passing
    # FALSE wrongly takes the "reuse old pANN" branch and errors out.
    seu <- doubletFinder(seu, PCs=1:10, pN=0.25, pK=best.pK, nExp=nExp_poi.adj, sct=FALSE)
    pann_col <- grep("^pANN", colnames(seu@meta.data), value=TRUE)
    seu@meta.data[[pann_col]]
  }
)

run_method <- function(method, counts, ds_name, seed) {
  st <- Sys.time()
  scores <- tryCatch(method_fns[[method]](counts, seed), error = function(e) {
    cat("    Failed on", method, "for", ds_name, "(seed=", seed, "):", conditionMessage(e), "\n")
    NULL
  })
  list(scores=scores, elapsed=as.numeric(Sys.time() - st, units="secs"))
}

new_results <- list()
new_percell <- list()
for (f in rds_files) {
  ds_name <- sub("\\.rds$", "", basename(f))

  obj <- readRDS(f)
  counts <- obj[[1]]
  truth_labels <- as.integer(tolower(obj[[2]]) == "doublet")
  cells <- colnames(counts)

  for (method in methods_to_run) {
    for (seed in SEEDS) {
      cat("Evaluating", ds_name, "-", method, "- seed", seed, "...\n")
      res <- run_method(method, counts, ds_name, seed)
      if (!is.null(res$scores)) {
        auroc <- as.numeric(pROC::auc(pROC::roc(truth_labels, res$scores, quiet=TRUE)))
        pr <- PRROC::pr.curve(scores.class0=res$scores, weights.class0=truth_labels)
        auprc <- pr$auc.integral

        new_results[[length(new_results) + 1]] <- data.frame(
          dataset=ds_name, method=method, seed=seed,
          AUPRC=auprc, AUROC=auroc, elapsed=res$elapsed
        )
        new_percell[[length(new_percell) + 1]] <- data.frame(
          dataset=ds_name, method=method, seed=seed,
          cell=cells, score=res$scores, truth=truth_labels
        )
      }
    }
  }
}

new_results_rows <- do.call(rbind, new_results)
new_percell_rows <- do.call(rbind, new_percell)

e <- rbind(kept_results, new_results_rows)
p <- rbind(kept_percell, new_percell_rows)

saveRDS(e, result_file)
saveRDS(p, percell_file)
write.csv(p, "benchmark_percell_scores_R.csv", row.names=FALSE)
cat("Saved R-package benchmark results to", result_file, "\n")
cat("Saved per-cell doublet scores to", percell_file, "and benchmark_percell_scores_R.csv\n")
