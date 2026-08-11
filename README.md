# scDblFinderPy paper

This is the data as well as the code to reproduce the figures from the scDblFinderPy paper, including the snakefile used to track the performance and resource usage of the different packages.

For the `scDblFinderPy` package itself, visit the package's [github page](https://github.com/ETHZ-INS/scDblFinderPy).

## Reproducing the Python benchmark

`benchmarking/run_python_benchmark.py` runs `scDblFinderPy` (in both
clustered and random mode) over each benchmark dataset, on 3 seeds by
default (`42, 0, 1` — pass `--seeds 1,2,3` to change them), and reports
AUPRC/AUROC/runtime per seed, matching the metrics used in the scDblFinder
paper's comparison figure. It also writes every cell's individual doublet
score to a separate CSV, which is what lets you later check whether
averaging the score across seeds (cell by cell) beats just averaging the
AUPRC of the 3 individual runs — see step 6 below.

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

Pass `--seeds` to change which seeds are run (default `42,0,1`):

```bash
python run_python_benchmark.py --seeds 1,2,3
```

For each `.h5ad` file in `benchmarking/datasets` and each seed, this
evaluates both clustered mode (`scDblFinder.Py.clusters`) and random mode
(`scDblFinder.Py.random`), and writes the resulting AUPRC/AUROC/runtime for
every dataset/method/seed row to
`benchmarking/python_benchmark_metrics_CPU.csv` (or `_GPU.csv` with
`--gpu`), plus every individual cell's doublet score to
`benchmarking/python_percell_scores_CPU.csv` (or `_GPU.csv`).

### 4. Rerun the other benchmarked methods

`benchmarking/plot_benchmark.R` (next step) compares `scDblFinderPy`
against Scrublet, Vaeda, `pyscdblfinder`, and a set of R-native methods, so
it needs their results on disk too. Each has its own script, and each also
runs 3 seeds by default (`--seeds` to override, same as step 3):

```bash
python run_scrublet_benchmark.py --seeds 42,0,1       # -> scrublet_benchmark_metrics.csv + scrublet_percell_scores.csv
/path/to/vaeda_env/bin/python run_vaeda_benchmark.py --seeds 42,0,1  # -> vaeda_benchmark_metrics.csv + vaeda_percell_scores.csv
python run_py_scdblfinder_benchmark.py --seeds 42,0,1 # -> py_scdblfinder_benchmark_metrics.csv + py_scdblfinder_percell_scores.csv
Rscript run_r_benchmark.R --all --seeds 42,0,1         # -> benchmark.results_R.rds + benchmark.percell_scores_R.rds (all R-native methods)
```

These need their own packages/environments (DoubletFinder, Seurat, scds,
the R `scDblFinder` package, pROC/PRROC, the `vaeda_env` conda env, and
`pyscdblfinder`) — see "Additional prerequisites" under the Snakemake
section below, which covers exactly the same scripts.

Every method here — including ones like DoubletFinder, bcds, and
`computeDoubletDensity` that don't expose an explicit seed argument — is
seeded via `set.seed()` (and, for DoubletFinder's Seurat steps, the
`seed.use`/`random.seed` arguments those functions default to fixed
constants for) so its 3 runs actually differ.

`run_r_benchmark.R` writes `benchmark.results_R.rds` (and
`benchmark.percell_scores_R.rds` for per-cell scores) and reads both back
on its next run to carry forward whatever it isn't asked to recompute, so
**the very first time you run it, pass `--all`** so every method has a
real value in there; after that, `Rscript run_r_benchmark.R` (no
arguments) reruns just `DoubletFinder` and keeps the rest as they were,
or pass a comma-separated list (e.g. `bcds,cxds`) to refresh a different
subset. See "Rerunning the other R-native methods" under the Snakemake
section below for more on this.

### 5. (Optional) Regenerate the comparison figures

`benchmarking/plot_benchmark.R` merges all of the above - both
`python_benchmark_metrics_CPU.csv` and `_GPU.csv` from step 3/4 (labelled
`scDblFinder.Py.clusters/.random (CPU)` and `(GPU)` respectively so they
appear as separate rows), `scrublet_benchmark_metrics.csv`,
`vaeda_benchmark_metrics.csv`, `py_scdblfinder_benchmark_metrics.csv`, plus
`benchmark.results_R.rds` - and produces two figures:

```bash
Rscript plot_benchmark.R
```

- `benchmark.pdf` — the main comparison figure, same layout as before, but
  now every method's AUPRC/AUROC/elapsed is first **averaged across its 3
  seeds** rather than plotting a single run.
- `benchmark_implementation_comparison.pdf` — a separate figure directly
  comparing the two independent Python reimplementations of scDblFinder
  (this repo's `scDblFinderPy`, both CPU and GPU, vs. the unrelated
  `pyscdblfinder` port, which only ever runs on CPU) seed by seed, as a
  scatter against the y=x parity line, faceted by clustered/random mode.

### 6. (Optional) Check whether averaging seeds actually helps

`benchmarking/analyze_averaging_gain.R` uses the per-cell score CSVs from
steps 3/4 to answer a different question: is the AUPRC computed from the
**cell-wise mean score across the 3 seeds** actually better than just the
**average of the 3 seeds' individual AUPRCs**? I.e., is there a real gain
from running a method multiple times and averaging, per (dataset, method):

```bash
Rscript analyze_averaging_gain.R
```

Produces `averaging_gain_summary.csv` (one row per dataset/method with both
numbers and their difference) and `benchmark_averaging_gain.pdf` (a
scatter of the same two numbers against the parity line, with a paired
Wilcoxon test printed to the console and in the plot subtitle).

## Automated pipeline (Snakemake)

`monitoring/Snakefile` automates steps 3-6 above for both a CPU and a GPU
run (including rerunning Scrublet, Vaeda, `pyscdblfinder`, and
DoubletFinder with their current versions, rather than relying solely on
the scDblFinder paper's original results), on 3 seeds per method by
default, tracks the resource usage (wall time, memory, CPU load, and — for
the GPU run — `nvidia-smi` utilization/memory) of each run, and produces
the main seed-averaged comparison figure, the scDblFinderPy-vs-pyscdblfinder
implementation-comparison figure, and the seed-averaging-gain analysis.

Steps 1-2 above (installing `scDblFinderPy` and fetching/converting the
datasets) are still prerequisites; Snakemake orchestrates the rest.

### 1. Install Snakemake

```bash
pip install snakemake
```

### 2. Additional prerequisites

On top of `scDblFinderPy` itself, a few more packages/environments are
needed to actually (re)compute each method rather than just carry its
previous value forward. Every method below runs from
`benchmarking/run_r_benchmark.R` unless noted otherwise, and that script
unconditionally loads `SingleCellExperiment`, `scds`, `scDblFinder`,
`Seurat`, and `DoubletFinder`, plus `pROC`/`PRROC` for scoring — so all of
these are required regardless of which subset of methods you actually
rerun via `rerun_methods` (see below).

- **DoubletFinder** (R, not on CRAN/Bioconductor — reruns the
  `DoubletFinder` method):
  ```bash
  Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
  ```
  This also pulls in **Seurat** (CRAN) as a dependency, which
  `run_r_benchmark.R` uses directly to build the object DoubletFinder
  expects.
- **scds** (R/Bioconductor — provides `bcds`/`cxds`/`hybrid`):
  ```bash
  Rscript -e 'BiocManager::install("scds", update=FALSE, ask=FALSE)'
  ```
- **scDblFinder** (R/Bioconductor — provides both the `scDblFinder.clusters`
  /`scDblFinder.random` methods *and* `computeDoubletDensity`, which now
  lives in this package rather than `scran`):
  ```bash
  Rscript -e 'BiocManager::install("scDblFinder", update=FALSE, ask=FALSE)'
  ```
- **pROC** and **PRROC** (CRAN — compute AUROC/AUPRC for every R-native
  method's scores):
  ```bash
  Rscript -e 'install.packages(c("pROC", "PRROC"))'
  ```
- **Scrublet** needs no extra installation: `benchmarking/run_scrublet_benchmark.py`
  uses `scanpy`'s built-in `sc.pp.scrublet`, and `scanpy` is already a
  dependency of `scDblFinderPy` from step 1 above.
- **pyscdblfinder** ([omicverse/py-scDblFinder](https://github.com/omicverse/py-scDblFinder),
  a separate pure-Python/xgboost port of the R `scDblFinder` algorithm,
  distinct from this repo's own `scDblFinderPy`) — install it into the
  same virtual environment used for `run_python_benchmark.py`:
  ```bash
  pip install pyscdblfinder
  ```
  Its pinned dependencies (`numpy`, `scipy`, `pandas`, `anndata`,
  `scikit-learn`, `xgboost`) are loose lower bounds, so this coexists with
  `scDblFinderPy`'s own pins from step 1 without needing a separate env.
  `benchmarking/run_py_scdblfinder_benchmark.py` uses its `ScDblFinder`
  class in the default (non-cluster-aware) mode, matching the package's
  own quick-start example.
- **Vaeda** (a Python/TensorFlow tool despite the "R packages" framing of
  the rest of this benchmark). It pins very old, specific dependencies
  (Python 3.8, TensorFlow 2.8) that conflict with the environment used for
  `run_python_benchmark.py`/`run_scrublet_benchmark.py`, so it needs its
  own conda environment:
  ```bash
  conda create -n vaeda_env python=3.8 -c conda-forge --override-channels -y
  conda run -n vaeda_env pip install tensorflow==2.8.0 tensorflow-probability==0.16.0 \
      'scanpy[leiden]'==1.8.0 typing-extensions==3.7.4 absl-py==0.10 six==1.15.0 \
      wrapt==1.12.1 xlrd==1.2.0 protobuf==3.19.6 matplotlib==3.5.3 pandas==1.3.5
  conda run -n vaeda_env pip install -i https://test.pypi.org/simple/ vaeda==0.0.30
  ```
  `monitoring/Snakefile`'s `run_vaeda_benchmark` rule invokes this env by
  the absolute path `/home/ahiropedi/.conda/envs/vaeda_env/bin/python`
  rather than `conda run -n vaeda_env` (which was unreliable on the
  machine this was developed on — it silently fell back to a different
  env's `python`). **If your `vaeda_env` lives elsewhere, update that path
  in the `run_vaeda_benchmark` rule** before running the pipeline.

### 3. Run it

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
snakemake -n                                    # dry run: show what would execute
snakemake --cores 1 benchmarking/benchmark.pdf   # build just the main figure
snakemake --cores 1 --config seeds=1,2,3         # use different seeds (default: 42,0,1)
```

The GPU rule requires the same conda/RAPIDS setup as `python
run_python_benchmark.py --gpu` above, plus `nvidia-smi` on `PATH`.

#### Rerunning the other R-native methods

`benchmarking/benchmark.results_R.rds` is self-perpetuating: each run of
`run_r_benchmark.R` reads whatever it last wrote and carries forward any
method it isn't asked to recompute this time (`Chord` is always dropped
as unmaintained, and `Scrublet`/`Vaeda` always come from their own Python
reruns rather than this file). **The first time you build it, use
`rerun_methods=all`** so every method gets a real value; after that, the
default reruns just `DoubletFinder` and keeps the rest as they were. Pass
`rerun_methods` as a comma-separated list or `all` to refresh more:

```bash
snakemake --cores 1 --config rerun_methods=all           # first run: compute every R-native method
snakemake --cores 1 --config rerun_methods=bcds,cxds     # later: refresh just a subset
```

### What it produces

- `benchmarking/python_benchmark_metrics_CPU.csv` / `_GPU.csv` — scDblFinderPy
  (clustered + random mode) results across the 3 seeds, same as the manual
  steps.
- `benchmarking/python_percell_scores_CPU.csv` / `_GPU.csv` — every cell's
  individual doublet score for each of the 3 seeds.
- `benchmarking/scrublet_benchmark_metrics.csv` /
  `scrublet_percell_scores.csv` — Scrublet, rerun fresh via
  `benchmarking/run_scrublet_benchmark.py`.
- `benchmarking/vaeda_benchmark_metrics.csv` / `vaeda_percell_scores.csv` —
  Vaeda, rerun fresh via `benchmarking/run_vaeda_benchmark.py` in the
  `vaeda_env` conda env.
- `benchmarking/py_scdblfinder_benchmark_metrics.csv` /
  `py_scdblfinder_percell_scores.csv` — `pyscdblfinder`, rerun fresh via
  `benchmarking/run_py_scdblfinder_benchmark.py`.
- `benchmarking/benchmark.results_R.rds` /
  `benchmarking/benchmark.percell_scores_R.rds` (+
  `benchmark_percell_scores_R.csv` export) — DoubletFinder rerun fresh
  (plus any other methods requested via `rerun_methods`) via
  `benchmarking/run_r_benchmark.R`, merged with whatever this same file
  already had for the methods that weren't rerun this time (see "Rerunning
  the other R-native methods" above — the very first run needs
  `rerun_methods=all`).
- `benchmarking/benchmark.pdf` — the main comparison figure, combining all
  of the above via `plot_benchmark.R`, with every method's AUPRC/AUROC/
  elapsed averaged across its 3 seeds first (Python CPU and GPU results
  shown as separate rows alongside every other method).
- `benchmarking/benchmark_implementation_comparison.pdf` — per-seed AUPRC
  scatter directly comparing `scDblFinderPy` (CPU and GPU) against
  `pyscdblfinder`, the other independent Python reimplementation of
  scDblFinder.
- `benchmarking/averaging_gain_summary.csv` and
  `benchmarking/benchmark_averaging_gain.pdf` — per (dataset, method)
  comparison of the AUPRC from cell-wise-averaging the score across the 3
  seeds vs. the average of the 3 seeds' individual AUPRCs, via
  `benchmarking/analyze_averaging_gain.R`.
- `monitoring/benchmarks/run_benchmark_CPU.tsv` / `_GPU.tsv` /
  `run_scrublet_benchmark.tsv` / `run_vaeda_benchmark.tsv` /
  `run_py_scdblfinder_benchmark.tsv` / `run_r_benchmark.tsv` — Snakemake's
  built-in per-run profile (wall time, max RSS, mean CPU load, I/O) for
  each of the above reruns.
- `monitoring/benchmarks/gpu_usage_GPU.csv` — GPU utilization % and memory
  used, sampled once a second for the duration of the GPU run.

The `monitoring/benchmarks/` outputs are run-machine-specific and are
gitignored rather than checked in.
