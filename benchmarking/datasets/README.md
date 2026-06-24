Populate this folder with:

```
curl -S https://zenodo.org/record/4562782/files/real_datasets.zip?download=1 > real_datasets.zip
unzip real_datasets.zip
```

Then convert the `.rds` files to `.h5ad` (requires R packages `SingleCellExperiment` and `zellkonverter`):

```r
BiocManager::install(c("SingleCellExperiment", "zellkonverter"))
```

```
Rscript convert_rds_to_h5ad.R
```