#!/usr/bin/env bash
# Run by Snakemake (--use-conda) right after creating the env from
# r_benchmark.yaml: installs the exact R package versions recorded in
# r_benchmark.renv.lock (compiled from source against this env's toolchain,
# so the first creation takes a while).
#
# Snakemake hashes this script's contents, not the lockfile's, to decide
# whether the env needs rebuilding - so the lockfile's checksum is embedded
# below, and must be updated whenever the lockfile is regenerated (the check
# fails loudly until it is).
set -euo pipefail

LOCK_SHA256=23dac169b498048f09424c44b6a321a3433c2f4e27c77adc8d8e0b3f6bc0a448
RENV_VERSION=1.3.1
CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-10-23

# Snakemake runs this from the workflow's workdir (the repo root).
lockfile=monitoring/envs/r_benchmark.renv.lock
echo "${LOCK_SHA256}  ${lockfile}" | sha256sum --check --quiet

# R invokes the compilers by bare name, so make sure this env's (which match
# the gcc r-base was built with) come first, whatever else is on PATH.
export PATH="$CONDA_PREFIX/bin:$PATH"

lib=$("$CONDA_PREFIX/bin/Rscript" -e 'cat(.libPaths()[1])')
# The renv cache would share compiled packages between envs, but they link
# against whichever env built them - so build everything inside this env.
export RENV_CONFIG_CACHE_ENABLED=FALSE
export MAKEFLAGS="-j${R_BUILD_JOBS:-8}"

"$CONDA_PREFIX/bin/Rscript" -e "
  options(repos=c(CRAN='${CRAN_SNAPSHOT}'), Ncpus=${R_BUILD_JOBS:-8})
  install.packages('remotes', lib='${lib}')
  remotes::install_version('renv', version='${RENV_VERSION}', lib='${lib}')
  renv::restore(lockfile='${lockfile}', library='${lib}', prompt=FALSE, clean=FALSE)
"
