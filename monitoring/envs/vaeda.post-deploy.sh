#!/usr/bin/env bash
# Run by Snakemake (--use-conda) right after creating the env from vaeda.yaml.
set -euo pipefail
"$CONDA_PREFIX/bin/python" -m pip install --no-deps --index-url https://test.pypi.org/simple/ vaeda==0.0.30
