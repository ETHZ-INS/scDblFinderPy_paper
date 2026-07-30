import os

# vaeda's pinned TensorFlow 2.8 crashes on this machine's GPU/thread setup
# (pthread_create failures against the current driver) unless GPUs are
# hidden and its thread pools are capped; must be set before importing tf.
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("TF_NUM_INTEROP_THREADS", "2")
os.environ.setdefault("TF_NUM_INTRAOP_THREADS", "2")

import scanpy as sc
import pandas as pd
import time
import glob
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

import vaeda


def eval_vaeda(adata, ds_name):
    st = time.time()
    try:
        res = vaeda.vaeda(adata.copy(), seed=42, verbose=0)
    except Exception as e:
        print(f"    Failed on Vaeda for {ds_name}: {e}")
        return None

    elapsed = time.time() - st

    truth_labels = (adata.obs['type'] == 'doublet').astype(int)
    scores = res.obs['vaeda_scores']

    auroc = roc_auc_score(truth_labels, scores)
    precision, recall, _ = precision_recall_curve(truth_labels, scores)
    auprc = auc(recall, precision)

    return {
        "dataset": ds_name,
        "method": "Vaeda",
        "AUPRC": auprc,
        "AUROC": auroc,
        "elapsed": elapsed
    }


def main():
    datasets = glob.glob(os.path.join(os.path.dirname(__file__), "datasets", "*.h5ad"))
    all_results = []

    print("Running Vaeda benchmark...")
    for ds in datasets:
        ds_name = os.path.basename(ds).replace(".h5ad", "")
        print(f"Evaluating {ds_name}...")

        adata = sc.read_h5ad(ds)
        if 'type' in adata.obs:
            adata.obs['type'] = adata.obs['type'].str.lower()

        res = eval_vaeda(adata, ds_name)
        if res:
            all_results.append(res)

    df = pd.DataFrame(all_results)
    out_path = os.path.join(os.path.dirname(__file__), "vaeda_benchmark_metrics.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved Vaeda benchmark scores across all datasets to {out_path}")


if __name__ == "__main__":
    main()
