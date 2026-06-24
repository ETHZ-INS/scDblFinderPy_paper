library(SingleCellExperiment)
library(zellkonverter)

rds_files <- list.files(".", pattern = "\\.rds$", full.names = TRUE)

if (length(rds_files) == 0) {
    stop("No .rds files found. Run the curl/unzip commands first.")
}

for (f in rds_files) {
    cat("Converting", f, "...\n")
    sce <- readRDS(f)
    out <- sub("\\.rds$", ".h5ad", f)
    writeH5AD(sce, out)
    cat("  -> Written to", out, "\n")
}
