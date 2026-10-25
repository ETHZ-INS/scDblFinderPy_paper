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

The datasets ship as `.rds`; they need converting to `.h5ad` (which
`run_python_benchmark.py` reads) with the provided `convert_rds_to_h5ad.R`.
The Snakemake pipeline below does this automatically (rule
`convert_datasets`, in the pinned R environment) whenever the `.h5ad` files
are missing. To do it by hand instead, run the script with the `Rscript` of
the pinned R environment (see "Pinned environments" below):

```bash
/path/to/r_benchmark_env/bin/Rscript convert_rds_to_h5ad.R
```

This conversion step can take a while given the number of datasets.

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
/path/to/scrublet_pyscdbl_env/bin/python run_scrublet_benchmark.py --seeds 42,0,1        # -> scrublet_benchmark_metrics.csv + scrublet_percell_scores.csv
/path/to/vaeda_env/bin/python run_vaeda_benchmark.py --seeds 42,0,1                      # -> vaeda_benchmark_metrics.csv + vaeda_percell_scores.csv
/path/to/scrublet_pyscdbl_env/bin/python run_py_scdblfinder_benchmark.py --seeds 42,0,1  # -> py_scdblfinder_benchmark_metrics.csv + py_scdblfinder_percell_scores.csv
/path/to/r_benchmark_env/bin/Rscript run_r_benchmark.R --all --seeds 42,0,1             # -> benchmark.results_R.rds + benchmark.percell_scores_R.rds (all R-native methods)
```

Each of these runs in its own pinned environment, separate from the
`scDblFinderPy` one — see "Pinned environments" under the Snakemake
section below for what they contain and how to build them (Snakemake
builds them automatically).

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
/path/to/r_benchmark_env/bin/Rscript plot_benchmark.R
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
/path/to/r_benchmark_env/bin/Rscript analyze_averaging_gain.R
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

Installing `scDblFinderPy` (step 1) and downloading the datasets (the
`curl`/`unzip` part of step 2) are still prerequisites; Snakemake
orchestrates the rest, including the `.rds` → `.h5ad` conversion.

### 1. Install Snakemake and conda

You need Snakemake ≥ 8 and conda (Miniconda, Miniforge, …):

```bash
pip install snakemake
```

Nothing else needs installing by hand: every method other than
`scDblFinderPy` runs in a pinned environment that Snakemake builds itself
on the first run with `--use-conda` (see below).

### 2. Pinned environments

Every rule except the two `scDblFinderPy` ones (`run_benchmark_cpu`/`_gpu`)
declares a conda environment from `monitoring/envs/`, with every package
version pinned so that anyone rerunning the pipeline gets exactly the
software the published results were produced with:

| Environment | Used by | What's pinned |
| --- | --- | --- |
| `r_benchmark.yaml` + `r_benchmark.renv.lock` | `run_r_benchmark.R`, `plot_benchmark.R`, `analyze_averaging_gain.R`, `convert_rds_to_h5ad.R` | R 4.6.1 and its compilers via conda; every R package (CRAN from a dated Posit Package Manager snapshot, Bioconductor 3.23, and `scDblFinder` / `DoubletFinder` at specific GitHub commits) via the renv lockfile |
| `scrublet_pyscdbl.yaml` | `run_scrublet_benchmark.py`, `run_py_scdblfinder_benchmark.py` | Python 3.13, `scanpy[scrublet]` 1.12.4 (Scrublet is `sc.pp.scrublet`), `pyscdblfinder` 0.3.0, plus every dependency |
| `vaeda.yaml` | `run_vaeda_benchmark.py` | Python 3.8, TensorFlow 2.8, plus every dependency; `vaeda` 0.0.30 itself comes from TestPyPI via `vaeda.post-deploy.sh` (installed with `--no-deps`, so nothing else is ever pulled from TestPyPI) |

Notes:

- **First run takes a while.** conda-forge doesn't ship R 4.6 builds of
  the Bioconductor stack yet, so `r_benchmark.post-deploy.sh` compiles every
  R package in the lockfile from source when Snakemake first creates that
  env. Set `R_BUILD_JOBS` (default 8) to change how many compile jobs run in
  parallel. The envs are stored under `.snakemake/conda/` and reused on
  every later run.
- `scDblFinder` is taken from its GitHub `devel` branch (pinned commit)
  rather than the Bioconductor release, so the R reference implementation
  matches the algorithm `scDblFinderPy` reimplements (e.g. the decaying
  doublet-rate weighting).
- To deliberately move to newer R package versions, edit the constants at
  the top of `monitoring/envs/r_benchmark_make_lock.R`, run it as described
  in its header to regenerate the lockfile, and update `LOCK_SHA256` (and
  `RENV_VERSION` if renv changed) in `r_benchmark.post-deploy.sh` — the
  post-deploy script refuses to run until the checksum matches, which is
  also what makes Snakemake rebuild the env.
- To use these envs outside Snakemake (the manual steps above), build them
  directly, e.g. from the repo root:
  ```bash
  conda env create -p envs/r_benchmark -f monitoring/envs/r_benchmark.yaml
  conda run -p envs/r_benchmark bash monitoring/envs/r_benchmark.post-deploy.sh
  ```
  (likewise for `vaeda`; `scrublet_pyscdbl` has no post-deploy step).
- If conda fails with `CondaToSNonInteractiveError`, your conda config
  still lists Anaconda's `defaults` channel, whose Terms of Service it
  wants accepted even though these envs only use conda-forge. Either accept
  them (`conda tos accept ...`, as the error suggests) or point `defaults`
  at conda-forge for the run:
  ```bash
  export CONDA_DEFAULT_CHANNELS=https://conda.anaconda.org/conda-forge
  ```

The two `scDblFinderPy` rules instead use whichever environment is active
when you invoke `snakemake`, so **activate the env from step 1 first** (for
the GPU rule, the conda/RAPIDS setup from the package's README, plus
`nvidia-smi` on `PATH`).

### 3. Run it

From the repo root, with the `scDblFinderPy` env active:

```bash
snakemake --cores 1 --use-conda -s monitoring/Snakefile
```

or, equivalently, from `monitoring/`:

```bash
cd monitoring
snakemake --cores 1 --use-conda
```

`--use-conda` is required: the Snakefile refuses to run without it, rather
than silently running every method in whatever env happens to be active.

Useful variations:

```bash
snakemake -n --use-conda                                    # dry run: show what would execute
snakemake --cores 1 --use-conda benchmarking/benchmark.pdf   # build just the main figure
snakemake --cores 1 --use-conda --config seeds=1,2,3         # use different seeds (default: 42,0,1)
snakemake --cores 1 --use-conda --conda-create-envs-only     # just build the pinned envs
```

#### Forcing a full rerun

Snakemake only reruns a rule when its inputs, code, parameters or conda env
changed. It does **not** see changes to the `scDblFinderPy` package itself
(it lives outside this repo), and existing outputs are otherwise left as
they are. To regenerate every result from scratch (without re-converting
the datasets):

```bash
snakemake --cores 1 --use-conda -s monitoring/Snakefile --rerun-incomplete \
  --forcerun run_benchmark_cpu run_benchmark_gpu run_scrublet_benchmark \
             run_py_scdblfinder_benchmark run_vaeda_benchmark run_r_benchmark \
  --config rerun_methods=all
```

The figures and the averaging-gain analysis rerun automatically after.
`--rerun-incomplete` makes Snakemake regenerate (rather than refuse to
start over) any outputs left behind by an earlier interrupted run; it's
harmless otherwise. The whole run takes many hours, so start it inside
`tmux`/`screen` (or with `nohup`) so it survives a dropped SSH session.

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
snakemake --cores 1 --use-conda --config rerun_methods=all           # first run: compute every R-native method
snakemake --cores 1 --use-conda --config rerun_methods=bcds,cxds     # later: refresh just a subset
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
  Vaeda, rerun fresh via `benchmarking/run_vaeda_benchmark.py` in its own
  pinned env.
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
