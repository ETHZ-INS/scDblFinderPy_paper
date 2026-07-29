# scDblFinderPy paper

This is the code & data behind the scDblFinderPy paper, as well as the snakefile used to track the performance and resource usage of the scDblFinderPy package.

For the `scDblFinderPy` package itself, visit the package's [github page](https://github.com/ETHZ-INS/scDblFinderPy).

## Reproducing the Python benchmark

`benchmarking/run_python_benchmark.py` runs `scDblFinderPy` (in both
clustered and random mode) over each benchmark dataset and reports
AUPRC/AUROC/runtime, matching the metrics used in the paper's comparison
figure.

### 1. Install scDblFinderPy

```bash
git clone https://github.com/ETHZ-INS/scDblFinderPy.git
cd scDblFinderPy
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install .
```

Requires Python 3.12+ (pinned dependencies such as `scanpy==1.12.1` and
`numpy==2.4.3` don't publish wheels for older versions). See the package's
[README](https://github.com/ETHZ-INS/scDblFinderPy#readme) for the optional
GPU setup (`use_gpu=True`), which needs a conda-built RAPIDS stack rather
than a plain pip install.

### 2. Get the benchmark datasets

From this repo's `benchmarking/datasets` directory:

```bash
cd benchmarking/datasets
curl -L -S "https://zenodo.org/record/4562782/files/real_datasets.zip?download=1" -o real_datasets.zip
unzip real_datasets.zip
```

The datasets ship as `.rds`; convert them to `.h5ad` (which
`run_python_benchmark.py` reads) with the provided R script. This requires
the R packages `SingleCellExperiment` and `zellkonverter`:

```bash
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("SingleCellExperiment", "zellkonverter"), update = FALSE, ask = FALSE)'
Rscript convert_rds_to_h5ad.R
```

This conversion step can take a while given the number of datasets and R
dependencies involved.

### 3. Run the benchmark

From the `benchmarking` directory, with the `scDblFinderPy` virtual
environment active:

```bash
cd ..   # back in benchmarking/
python run_python_benchmark.py
```

Pass `--gpu` to run with GPU acceleration (requires the conda/RAPIDS setup
from step 1):

```bash
python run_python_benchmark.py --gpu
```

For each `.h5ad` file in `benchmarking/datasets`, this evaluates both
clustered mode (`scDblFinder.Py.clusters`) and random mode
(`scDblFinder.Py.random`), and writes the resulting AUPRC/AUROC/runtime
for every dataset/method pair to `benchmarking/python_benchmark_metrics.csv`.

### 4. (Optional) Regenerate the comparison figure

`benchmarking/plot_benchmark.R` merges `python_benchmark_metrics.csv` with
the paper's original R-method results (`benchmark.results.rds`) to produce
the AUPRC comparison figure (`benchmark_AUPRC_fig.png`):

```bash
Rscript plot_benchmark.R
```

## Automated pipeline (Snakemake)

`monitoring/Snakefile` automates steps 3-4 above for both a CPU and a GPU
run, and additionally tracks the resource usage (wall time, memory, CPU
load, and — for the GPU run — `nvidia-smi` utilization/memory) of each run.

Steps 1-2 above (installing `scDblFinderPy` and fetching/converting the
datasets) are still prerequisites; Snakemake orchestrates the rest.

### 1. Install Snakemake

```bash
pip install snakemake
```

### 2. Run it

From the repo root:

```bash
snakemake --cores 1 -s monitoring/Snakefile
```

or, equivalently, from `monitoring/`:

```bash
cd monitoring
snakemake --cores 1
```

Useful variations:

```bash
snakemake -n                                             # dry run: show what would execute
snakemake --cores 1 benchmarking/benchmark_AUPRC_fig_GPU.png  # build just one target
```

The GPU rule requires the same conda/RAPIDS setup as `python
run_python_benchmark.py --gpu` above, plus `nvidia-smi` on `PATH`.

### What it produces

- `benchmarking/python_benchmark_metrics_CPU.csv` / `_GPU.csv` and
  `benchmarking/benchmark_AUPRC_fig_CPU.png` / `_GPU.png` — same outputs as
  the manual steps, just written straight to their CPU/GPU-suffixed names.
- `monitoring/benchmarks/run_benchmark_CPU.tsv` / `_GPU.tsv` — Snakemake's
  built-in per-run profile (wall time, max RSS, mean CPU load, I/O).
- `monitoring/benchmarks/gpu_usage_GPU.csv` — GPU utilization % and memory
  used, sampled once a second for the duration of the GPU run.

The `monitoring/benchmarks/` outputs are run-machine-specific and are
gitignored rather than checked in.