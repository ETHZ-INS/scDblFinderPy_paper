# scDblFinderPy paper

This is the code & data behind the scDblFinderPy paper.

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