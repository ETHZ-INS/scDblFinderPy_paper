library(ggplot2)
library(viridisLite)
library(reshape2)
library(patchwork)

cat("Loading R-package benchmark results...\n")
e <- readRDS("benchmark.results_R.rds")

cat("Loading Python package performance (CPU + GPU)...\n")
# Both CSVs use the same method names (scDblFinder.Py.clusters/.random), so
# tag each with its device before combining - otherwise the CPU and GPU rows
# for the same method would collide into a single row in the plot below.
py_cpu_res <- read.csv("python_benchmark_metrics_CPU.csv")
py_cpu_res$method <- paste0(py_cpu_res$method, " (CPU)")
py_gpu_res <- read.csv("python_benchmark_metrics_GPU.csv")
py_gpu_res$method <- paste0(py_gpu_res$method, " (GPU)")

scrublet_res <- read.csv("scrublet_benchmark_metrics.csv")
vaeda_res <- read.csv("vaeda_benchmark_metrics.csv")
py_scdblfinder_res <- read.csv("py_scdblfinder_benchmark_metrics.csv")
py_res <- rbind(py_cpu_res, py_gpu_res, scrublet_res, vaeda_res, py_scdblfinder_res)
# Rename columns if needed to match e
# e has: dataset, method, seed, AUPRC, AUROC, elapsed

for (i in 1:nrow(py_res)) {
  e <- rbind(e, c(py_res$dataset[i], py_res$method[i], py_res$seed[i],
                   py_res$AUPRC[i], py_res$AUROC[i], py_res$elapsed[i]))
}

# Ensure formats
e$seed <- as.integer(e$seed)
e$AUPRC <- as.numeric(e$AUPRC)
e$AUROC <- as.numeric(e$AUROC)
e$elapsed <- as.numeric(e$elapsed)

# Every method here is run on 3 seeds (see run_r_benchmark.R / run_python_
# benchmark.py etc.). Keep the full per-seed table around for the separate
# implementation-comparison figure below, but the main figure below
# summarises across seeds first, i.e. it plots the *average* AUPRC/AUROC
# per (dataset, method) - matching Pierre-Luc's request that the main
# figure show the seed-averaged AUPRC, not any single seed's value.
e_perseed <- e

cat("Averaging AUPRC/AUROC/elapsed across seeds for the main figure...\n")
e <- aggregate(cbind(AUPRC, AUROC, elapsed) ~ dataset + method, data=e, FUN=mean)

cat("Calculating metrics and ranks (as per paper script)...\n")
datmax <- sort(apply(reshape2::dcast(e, method~dataset, value.var="AUPRC")[,-1], 2, na.rm=TRUE, FUN=max))
metmax <- sort(apply(reshape2::dcast(e, dataset~method, value.var="AUPRC")[,-1], 2,
                     FUN=function(x){ mean(x, na.rm=TRUE) }))

met_levels <- names(metmax)

e$dataset <- factor(e$dataset, levels=names(datmax))
e$method <- factor(e$method, levels=met_levels)
levels(e$method) <- gsub("bcds","scds::bcds",levels(e$method))
levels(e$method) <- gsub("cxds","scds::cxds",levels(e$method))
levels(e$method) <- gsub("hybrid","scds::hybrid",levels(e$method))
levels(e$method) <- gsub("\\.Py","Py",levels(e$method))
levels(e$method) <- gsub("pyscDblFinder","pyscDblFinder.random", levels(e$method), fixed=TRUE)

getranks <- function(x){
  y <- rank(x)
  y[is.na(x)] <- NA
  y
}

tr <- reshape2::dcast(e, method~dataset, value.var="AUPRC")
row.names(tr) <- tr[,1]; tr <- tr[,-1]
tr2 <- apply(tr, 2, FUN=getranks)

e$AUPRC.rank <- apply(e[,1:2], 1, FUN=function(x) tr2[as.character(x[2]),as.character(x[1])])
e$point.colour <- viridisLite::viridis(100)[round(100*e$AUPRC)]
e$border.colour <- ifelse(e$AUPRC.rank>=(length(metmax)-0.5), "black", NA)
e$rounded <- round(e$AUPRC,2)
e$text.colour <- ifelse(e$rounded>=0.5,"black","white")
e$text.colour[e$AUPRC.rank==1] <- NA
e$text <- gsub("1\\.00","1.0",gsub("0\\.",".",sprintf("%.2f",e$AUPRC)))

p1 <- ggplot(e, aes(dataset, method)) +
  geom_point(aes(size=AUPRC, colour=AUPRC.rank)) +
  scale_size(range=c(4,10)) +
  scale_colour_viridis_c(breaks=c(1, length(unique(e$method))), labels=c("worst","best"), guide=guide_colorbar(barheight=unit(1, "cm"))) +
  geom_point(data=e[e$AUPRC.rank==length(unique(e$method)),], shape=21, colour="black", aes(size=AUPRC), stroke=1.1, show.legend=FALSE) +
  geom_text(aes(label=text), size=3, colour=ifelse(e$AUPRC.rank>=5, "black", "white")) +
  labs(colour="AUPRC rank", x="Dataset") + theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        axis.text.y=element_text(hjust=0.5, size=10.5, face="bold",
                                 colour=ifelse(levels(e$method) %in% scdbl.methods,"black","grey30")),
        axis.title.y=element_blank(), panel.grid=element_blank())

p2 <- ggplot(e, aes(elapsed, method)) + geom_col(fill="darkblue") + theme_minimal() +
  theme(axis.text.y=element_blank(), axis.ticks.x=element_line(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_x_reverse(breaks=scales::breaks_pretty(n=3)) +
  labs(x="Total time (seconds)", y="Method")

pp <- p2 + p1 + plot_layout(widths=c(1,5))

ggsave("benchmark.pdf", pp, width=11, height=6, dpi=300)
#ggsave("benchmark_AUPRC_fig.png", pp, width=12, height=6.5, dpi=300)
cat("Main (seed-averaged) plot generated successfully!\n")

# --- Implementation comparison figure -------------------------------------
# Directly compares the two independent Python reimplementations of
# scDblFinder (this repo's scDblFinderPy vs. the separate pyscDblFinder
# port, https://github.com/omicverse/py-scDblFinder) seed by seed, per
# Pierre-Luc's follow-up: "if the current results hold, we could then in a
# separate figure use the results of the individual seeds to directly
# compare the two implementations." pyscDblFinder has no clusters/random
# distinction (one score per dataset/seed) and, despite its 'auto' PCA/kNN
# backends, no actual working GPU path end-to-end (its xgboost classifier
# has an internal device= knob that's never wired up through the public
# API) - so it only ever runs on CPU and is compared as such against BOTH
# of scDblFinder.Py's CPU and GPU runs, one point per (dataset, mode,
# device, seed); x=y line marks parity between implementations.
cat("Building implementation-comparison figure (scDblFinderPy vs. pyscDblFinder, per seed)...\n")

pysc_impl <- e_perseed[e_perseed$method == "pyscDblFinder", c("dataset", "seed", "AUPRC")]
names(pysc_impl)[names(pysc_impl) == "AUPRC"] <- "AUPRC_pyscDblFinder"

py_impl <- e_perseed[e_perseed$method %in% c(
  "scDblFinder.Py.clusters (CPU)", "scDblFinder.Py.random (CPU)",
  "scDblFinder.Py.clusters (GPU)", "scDblFinder.Py.random (GPU)"
), ]
py_impl$mode <- ifelse(grepl("clusters", py_impl$method, fixed=TRUE), "clusters", "random")
py_impl$device <- ifelse(grepl("(CPU)", py_impl$method, fixed=TRUE), "CPU", "GPU")
py_impl <- py_impl[, c("dataset", "mode", "device", "seed", "AUPRC")]
names(py_impl)[names(py_impl) == "AUPRC"] <- "AUPRC_scDblFinderPy"

cmp <- merge(pysc_impl, py_impl, by=c("dataset", "seed"))

p3 <- ggplot(cmp, aes(AUPRC_pyscDblFinder, AUPRC_scDblFinderPy, colour=dataset, shape=device)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
  geom_point(size=2.5, alpha=0.85) +
  facet_wrap(~mode) +
  coord_equal(xlim=c(0,1), ylim=c(0,1)) +
  labs(x="AUPRC - pyscDblFinder (CPU)", y="AUPRC - scDblFinder.Py",
       title="Per-seed AUPRC: scDblFinder.Py vs. pyscDblFinder",
       colour="Dataset", shape="scDblFinder.Py device") +
  theme_minimal() +
  theme(legend.position="right")

ggsave("benchmark_implementation_comparison.pdf", p3, width=10, height=5, dpi=300)
cat("Implementation-comparison plot generated successfully!\n")
