# Regenerates r_benchmark.renv.lock - only needed when deliberately bumping
# the R package versions, not for running the pipeline (which restores the
# committed lockfile via r_benchmark.post-deploy.sh instead).
#
# Run with the Rscript of a fresh env built from r_benchmark.yaml, e.g.:
#   conda env create -p ./r_boot -f monitoring/envs/r_benchmark.yaml
#   PATH="$PWD/r_boot/bin:$PATH" ./r_boot/bin/Rscript monitoring/envs/r_benchmark_make_lock.R
#
# CRAN packages come from a dated Posit Package Manager snapshot and
# Bioconductor ones from a fixed release, so rerunning this script later
# resolves to the same versions unless the constants below are changed.

CRAN_SNAPSHOT <- "https://packagemanager.posit.co/cran/2026-10-23"
BIOC_VERSION <- "3.23"
# Latest commits at the time of the benchmark rerun; scDblFinder is taken
# from its devel branch rather than the Bioconductor release so that it
# matches the algorithm scDblFinderPy reimplements (e.g. dbr decay).
GITHUB_PKGS <- c(
  "plger/scDblFinder@1a0831bd096dad402af8c3b3e4d6e472b88ffab8",
  "chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d54b4b1cb4365639931bbb16f01e1cd"
)
# Every package the benchmarking/ R scripts load, plus BiocManager (which
# renv needs to restore Bioconductor packages).
TOP_LEVEL_PKGS <- c(
  "BiocManager", "remotes",
  "SingleCellExperiment", "scds", "scDblFinder", "Seurat", "DoubletFinder",
  "pROC", "PRROC", "zellkonverter",
  "ggplot2", "patchwork", "reshape2", "scales", "viridisLite"
)

lockfile <- file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value=TRUE))), "r_benchmark.renv.lock")
lib <- .libPaths()[1]

options(repos=c(CRAN=CRAN_SNAPSHOT), Ncpus=16)
install.packages(c("BiocManager", "remotes", "renv"), lib=lib)
options(repos=BiocManager::repositories(version=BIOC_VERSION))

cran_bioc <- setdiff(TOP_LEVEL_PKGS, c("BiocManager", "remotes", "scDblFinder", "DoubletFinder"))
BiocManager::install(cran_bioc, version=BIOC_VERSION, lib=lib,
                     update=FALSE, ask=FALSE)
# Installs any dependencies of the GitHub packages that the line above
# didn't already pull in (e.g. scDblFinder's own Bioconductor imports).
remotes::install_github(GITHUB_PKGS, lib=lib, upgrade="never",
                        repos=getOption("repos"))

# install.packages() only warns when a package fails to build, which would
# otherwise leave it (and everything depending on it) out of the lockfile.
missing <- setdiff(TOP_LEVEL_PKGS, rownames(installed.packages(lib.loc=lib)))
if (length(missing) > 0) {
  stop("Failed to install: ", paste(missing, collapse=", "))
}

project <- tempfile("r_benchmark_lock_")
dir.create(project)
renv::snapshot(project=project, library=lib, lockfile=lockfile,
               packages=TOP_LEVEL_PKGS, prompt=FALSE)
cat("Wrote", lockfile, "\n")
