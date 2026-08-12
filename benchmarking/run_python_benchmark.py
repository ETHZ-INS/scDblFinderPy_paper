import scanpy as sc
import pandas as pd
import time
import sys
import os
import gc
import glob
from sklearn.metrics import roc_auc_score, precision_recall_curve, auc


def parse_seeds(argv, default="42,0,1"):
    if '--seeds' in argv:
        raw = argv[argv.index('--seeds') + 1]
    else:
        raw = default
    return [int(s) for s in raw.split(',')]


from scDblFinderPy.scDblFinder import compute_doublet_score

def eval_mode(adata, mode_name, clusters_col, ds_name, seed, use_gpu=False):
    st = time.time()
    try:
        adata_res = compute_doublet_score(adata.copy(), n_iters=3, random_state=seed,
                                          clusters_col=clusters_col, use_gpu=use_gpu)
    except Exception as e:
        print(f"    Failed on {mode_name} for {ds_name} (seed={seed}): {e}")
        return None, None
    finally:
        # cupy's memory pools cache every block they ever allocate and never
        # release it back to the OS, so repeated GPU calls across the many
        # dataset/mode/seed combinations in this loop steadily inflate the
        # process's virtual memory until it hits the session's ulimit -v,
        # after which every subsequent dataset fails identically.
        if use_gpu:
            gc.collect()
            import cupy as cp
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()

    elapsed = time.time() - st

    truth_labels = (adata_res.obs['type'] == 'doublet').astype(int)
    scores = adata_res.obs['scDblFinder_score']

    auroc = roc_auc_score(truth_labels, scores)
    precision, recall, _ = precision_recall_curve(truth_labels, scores)
    auprc = auc(recall, precision)

    metrics_row = {
        "dataset": ds_name,
        "method": mode_name,
        "seed": seed,
        "AUPRC": auprc,
        "AUROC": auroc,
        "elapsed": elapsed
    }
    percell_rows = pd.DataFrame({
        "dataset": ds_name,
        "method": mode_name,
        "seed": seed,
        "cell": adata_res.obs_names,
        "score": scores.values,
        "truth": truth_labels.values,
    })
    return metrics_row, percell_rows

def main():
    use_gpu = '--gpu' in sys.argv
    seeds = parse_seeds(sys.argv)
    datasets = glob.glob(os.path.join(os.path.dirname(__file__), "datasets", "*.h5ad"))
    all_results = []
    all_percell = []

    print(f"Running benchmark (use_gpu={use_gpu}, seeds={seeds})...")
    for ds in datasets:
        ds_name = os.path.basename(ds).replace(".h5ad", "")
        print(f"Evaluating {ds_name}...")

        adata = sc.read_h5ad(ds)
        # some datasets string parsing
        if 'type' in adata.obs:
            adata.obs['type'] = adata.obs['type'].str.lower()

        for seed in seeds:
            print(f"  -> Clustered Mode (seed={seed})")
            res_clust, pc_clust = eval_mode(adata, "scDblFinder.Py.clusters", "clusters", ds_name,
                                 seed, use_gpu=use_gpu)
            if res_clust:
                all_results.append(res_clust)
                all_percell.append(pc_clust)

            print(f"  -> Random Mode (seed={seed})")
            res_rand, pc_rand = eval_mode(adata, "scDblFinder.Py.random", None, ds_name,
                                 seed, use_gpu=use_gpu)
            if res_rand:
                all_results.append(res_rand)
                all_percell.append(pc_rand)

    df = pd.DataFrame(all_results)
    suffix = "GPU" if use_gpu else "CPU"
    out_path = os.path.join(os.path.dirname(__file__), f"python_benchmark_metrics_{suffix}.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved python benchmark scores across all datasets to {out_path}")

    percell_df = pd.concat(all_percell, ignore_index=True)
    percell_path = os.path.join(os.path.dirname(__file__), f"python_percell_scores_{suffix}.csv")
    percell_df.to_csv(percell_path, index=False)
    print(f"Saved per-cell doublet scores to {percell_path}")

if __name__ == "__main__":
    main()
