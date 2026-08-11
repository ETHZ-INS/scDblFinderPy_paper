import scanpy as sc
import pandas as pd
import time
import sys
import os
import glob
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc


def parse_seeds(argv, default="42,0,1"):
    if '--seeds' in argv:
        raw = argv[argv.index('--seeds') + 1]
    else:
        raw = default
    return [int(s) for s in raw.split(',')]


def eval_scrublet(adata, ds_name, seed):
    st = time.time()
    try:
        work = adata.copy()
        sc.pp.scrublet(work, random_state=seed, verbose=False)
    except Exception as e:
        print(f"    Failed on Scrublet for {ds_name} (seed={seed}): {e}")
        return None, None

    elapsed = time.time() - st

    truth_labels = (adata.obs['type'] == 'doublet').astype(int)
    scores = work.obs['doublet_score']

    auroc = roc_auc_score(truth_labels, scores)
    precision, recall, _ = precision_recall_curve(truth_labels, scores)
    auprc = auc(recall, precision)

    metrics_row = {
        "dataset": ds_name,
        "method": "Scrublet",
        "seed": seed,
        "AUPRC": auprc,
        "AUROC": auroc,
        "elapsed": elapsed
    }
    percell_rows = pd.DataFrame({
        "dataset": ds_name,
        "method": "Scrublet",
        "seed": seed,
        "cell": work.obs_names,
        "score": scores.values,
        "truth": truth_labels.values,
    })
    return metrics_row, percell_rows


def main():
    seeds = parse_seeds(sys.argv)
    datasets = glob.glob(os.path.join(os.path.dirname(__file__), "datasets", "*.h5ad"))
    all_results = []
    all_percell = []

    print(f"Running Scrublet benchmark (seeds={seeds})...")
    for ds in datasets:
        ds_name = os.path.basename(ds).replace(".h5ad", "")
        print(f"Evaluating {ds_name}...")

        adata = sc.read_h5ad(ds)
        if 'type' in adata.obs:
            adata.obs['type'] = adata.obs['type'].str.lower()

        for seed in seeds:
            res, pc = eval_scrublet(adata, ds_name, seed)
            if res:
                all_results.append(res)
                all_percell.append(pc)

    df = pd.DataFrame(all_results)
    out_path = os.path.join(os.path.dirname(__file__), "scrublet_benchmark_metrics.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved Scrublet benchmark scores across all datasets to {out_path}")

    percell_df = pd.concat(all_percell, ignore_index=True)
    percell_path = os.path.join(os.path.dirname(__file__), "scrublet_percell_scores.csv")
    percell_df.to_csv(percell_path, index=False)
    print(f"Saved per-cell doublet scores to {percell_path}")


if __name__ == "__main__":
    main()
