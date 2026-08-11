suppressMessages({
  library(ggplot2)
  library(reshape2)
})

# Addresses Pierre-Luc's question: "is the AUPRC calculated on the cell-wise
# average score across the 3 seeds significantly better than the average
# AUPRC of the 3 seeds separately?" - i.e. is there a real gain from running
# a method multiple times (with different seeds) and averaging the scores,
# vs. just running it once. Uses the per-cell score files each benchmark
# script now writes (one row per dataset/method/seed/cell).

cat("Loading per-cell doublet scores...\n")

r_pc <- read.csv("benchmark_percell_scores_R.csv")

py_cpu_pc <- read.csv("python_percell_scores_CPU.csv")
py_cpu_pc$method <- paste0(py_cpu_pc$method, " (CPU)")
py_gpu_pc <- read.csv("python_percell_scores_GPU.csv")
py_gpu_pc$method <- paste0(py_gpu_pc$method, " (GPU)")

scrublet_pc <- read.csv("scrublet_percell_scores.csv")
vaeda_pc <- read.csv("vaeda_percell_scores.csv")
py_scdblfinder_pc <- read.csv("py_scdblfinder_percell_scores.csv")

pc <- rbind(r_pc, py_cpu_pc, py_gpu_pc, scrublet_pc, vaeda_pc, py_scdblfinder_pc)
pc$seed <- as.integer(pc$seed)
pc$score <- as.numeric(pc$score)
pc$truth <- as.integer(pc$truth)

auprc <- function(scores, truth) {
  pr <- PRROC::pr.curve(scores.class0=scores, weights.class0=truth)
  pr$auc.integral
}

cat("Computing per-seed AUPRC and cell-wise-averaged-score AUPRC for each (dataset, method)...\n")

results <- list()
for (key in split(pc, list(pc$dataset, pc$method), drop=TRUE)) {
  if (length(unique(key$seed)) < 2) next  # nothing to average
  ds_name <- key$dataset[1]
  method <- key$method[1]

  # Mean of the per-seed AUPRCs (i.e. run once, repeat, average the metric).
  per_seed_auprc <- sapply(split(key, key$seed), function(s) auprc(s$score, s$truth))
  avg_of_seed_auprcs <- mean(per_seed_auprc)

  # AUPRC of the cell-wise mean score across seeds (i.e. run N times,
  # average the scores per cell first, then compute one AUPRC).
  wide <- reshape2::dcast(key, cell + truth ~ seed, value.var="score")
  score_cols <- setdiff(colnames(wide), c("cell", "truth"))
  mean_score <- rowMeans(wide[, score_cols, drop=FALSE], na.rm=TRUE)
  auprc_of_avg_score <- auprc(mean_score, wide$truth)

  results[[length(results) + 1]] <- data.frame(
    dataset=ds_name, method=method, n_seeds=length(unique(key$seed)),
    avg_of_seed_AUPRCs=avg_of_seed_auprcs,
    AUPRC_of_avg_score=auprc_of_avg_score,
    gain=auprc_of_avg_score - avg_of_seed_auprcs
  )
}
summary_df <- do.call(rbind, results)

out_csv <- "averaging_gain_summary.csv"
write.csv(summary_df, out_csv, row.names=FALSE)
cat("Saved per-(dataset, method) comparison to", out_csv, "\n")

test <- wilcox.test(summary_df$AUPRC_of_avg_score, summary_df$avg_of_seed_AUPRCs, paired=TRUE)
cat(sprintf(
  "\nAcross %d (dataset, method) combinations: mean gain from cell-wise score averaging = %.4f AUPRC (paired Wilcoxon p = %.4g)\n",
  nrow(summary_df), mean(summary_df$gain), test$p.value
))

p <- ggplot(summary_df, aes(avg_of_seed_AUPRCs, AUPRC_of_avg_score)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey50") +
  geom_point(aes(colour=method), size=2, alpha=0.8) +
  coord_equal() +
  labs(x="Mean of the 3 seeds' individual AUPRCs",
       y="AUPRC of the 3-seed cell-wise-averaged score",
       title="Does averaging doublet scores across seeds improve AUPRC?",
       subtitle=sprintf("Paired Wilcoxon p = %.4g, mean gain = %.4f", test$p.value, mean(summary_df$gain)),
       colour="Method") +
  theme_minimal() +
  theme(legend.position="right")

ggsave("benchmark_averaging_gain.pdf", p, width=8, height=7, dpi=300)
cat("Saved plot to benchmark_averaging_gain.pdf\n")
