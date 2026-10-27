import os

# vaeda's pinned TensorFlow 2.8 crashes on this machine's GPU/thread setup
# (pthread_create failures against the current driver) unless GPUs are
# hidden and its thread pools are capped; must be set before importing tf.
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("TF_NUM_INTEROP_THREADS", "1")
os.environ.setdefault("TF_NUM_INTRAOP_THREADS", "1")

# Even with single-threaded pools, TensorFlow runs ops on a worker thread
# alongside the Python main thread, and the two overlap (~110% CPU), so pin
# the process (and every thread it starts) to one core to match the other
# methods' single-core benchmark.
os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})

import scanpy as sc
import pandas as pd
import time
import sys
import glob
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

import vaeda


def parse_seeds(argv, default="42,0,1"):
    if '--seeds' in argv:
        raw = argv[argv.index('--seeds') + 1]
    else:
        raw = default
    return [int(s) for s in raw.split(',')]


def eval_vaeda(adata, ds_name, seed):
    st = time.time()
    try:
        res = vaeda.vaeda(adata.copy(), seed=seed, verbose=0)
    except Exception as e:
        print(f"    Failed on Vaeda for {ds_name} (seed={seed}): {e}")
        return None, None

    elapsed = time.time() - st

    truth_labels = (adata.obs['type'] == 'doublet').astype(int)
    scores = res.obs['vaeda_scores']

    auroc = roc_auc_score(truth_labels, scores)
    precision, recall, _ = precision_recall_curve(truth_labels, scores)
    auprc = auc(recall, precision)

    metrics_row = {
        "dataset": ds_name,
        "method": "Vaeda",
        "seed": seed,
        "AUPRC": auprc,
        "AUROC": auroc,
        "elapsed": elapsed
    }
    percell_rows = pd.DataFrame({
        "dataset": ds_name,
        "method": "Vaeda",
        "seed": seed,
        "cell": res.obs_names,
        "score": scores.values,
        "truth": truth_labels.values,
    })
    return metrics_row, percell_rows


def main():
    seeds = parse_seeds(sys.argv)
    datasets = glob.glob(os.path.join(os.path.dirname(__file__), "datasets", "*.h5ad"))
    all_results = []
    all_percell = []

    print(f"Running Vaeda benchmark (seeds={seeds})...")
    for ds in datasets:
        ds_name = os.path.basename(ds).replace(".h5ad", "")
        print(f"Evaluating {ds_name}...")

        adata = sc.read_h5ad(ds)
        if 'type' in adata.obs:
            adata.obs['type'] = adata.obs['type'].str.lower()

        for seed in seeds:
            res, pc = eval_vaeda(adata, ds_name, seed)
            if res:
                all_results.append(res)
                all_percell.append(pc)

    df = pd.DataFrame(all_results)
    out_path = os.path.join(os.path.dirname(__file__), "vaeda_benchmark_metrics.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved Vaeda benchmark scores across all datasets to {out_path}")

    percell_df = pd.concat(all_percell, ignore_index=True)
    percell_path = os.path.join(os.path.dirname(__file__), "vaeda_percell_scores.csv")
    percell_df.to_csv(percell_path, index=False)
    print(f"Saved per-cell doublet scores to {percell_path}")


if __name__ == "__main__":
    main()
