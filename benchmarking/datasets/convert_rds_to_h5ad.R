library(SingleCellExperiment)
library(zellkonverter)

rds_files <- list.files(".", pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
    stop("No .rds files found. Run the curl/unzip commands first.")
}

for (f in rds_files) {
    cat("Converting", f, "...\n")
    obj <- readRDS(f)
    if (is(obj, "SingleCellExperiment")) {
        sce <- obj
    } else {
        # Zenodo real_datasets.zip stores each dataset as an unnamed
        # list(counts_matrix, truth_labels) rather than an SCE.
        counts <- obj[[1]]
        truth <- obj[[2]]
        sce <- SingleCellExperiment(
            assays = list(counts = counts),
            colData = DataFrame(type = truth)
        )
    }
    out <- sub("\\.rds$", ".h5ad", f)
    writeH5AD(sce, out)
    cat("  -> Written to", out, "\n")
}
