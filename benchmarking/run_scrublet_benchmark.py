import scanpy as sc
import pandas as pd
import time
import os
import glob
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc


def eval_scrublet(adata, ds_name):
    st = time.time()
    try:
        work = adata.copy()
        sc.pp.scrublet(work, random_state=42, verbose=False)
    except Exception as e:
        print(f"    Failed on Scrublet for {ds_name}: {e}")
        return None

    elapsed = time.time() - st

    truth_labels = (adata.obs['type'] == 'doublet').astype(int)
    scores = work.obs['doublet_score']

    auroc = roc_auc_score(truth_labels, scores)
    precision, recall, _ = precision_recall_curve(truth_labels, scores)
    auprc = auc(recall, precision)

    return {
        "dataset": ds_name,
        "method": "Scrublet",
        "AUPRC": auprc,
        "AUROC": auroc,
        "elapsed": elapsed
    }


def main():
    datasets = glob.glob(os.path.join(os.path.dirname(__file__), "datasets", "*.h5ad"))
    all_results = []

    print("Running Scrublet benchmark...")
    for ds in datasets:
        ds_name = os.path.basename(ds).replace(".h5ad", "")
        print(f"Evaluating {ds_name}...")

        adata = sc.read_h5ad(ds)
        if 'type' in adata.obs:
            adata.obs['type'] = adata.obs['type'].str.lower()

        res = eval_scrublet(adata, ds_name)
        if res:
            all_results.append(res)

    df = pd.DataFrame(all_results)
    out_path = os.path.join(os.path.dirname(__file__), "scrublet_benchmark_metrics.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved Scrublet benchmark scores across all datasets to {out_path}")


if __name__ == "__main__":
    main()
