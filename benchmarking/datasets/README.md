Populate this folder with:

```
curl -L -S "https://zenodo.org/record/4562782/files/real_datasets.zip?download=1" -o real_datasets.zip

unzip real_datasets.zip
```

Then convert the `.rds` files to `.h5ad` (requires R packages `SingleCellExperiment` and `zellkonverter`).:

```
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("SingleCellExperiment", "zellkonverter"), update = FALSE, ask = FALSE)'

Rscript convert_rds_to_h5ad.R
```

Note: the conversion process may be quite time consuming due to the large number of datasets and the need to set up all the dependencies correctly.
