rm(list=ls())
library(data.table)
library(ggplot2)
source("paper/functions.R")

validation_path <- "data/validation/"
fig_path <- "paper/figures/simul_supp/"
if (!dir.exists(fig_path)) {
  dir.create(fig_path, recursive = T)
}

method_levels <- c("DESeq2", "EdgeR", "NBSR")
cols <- RColorBrewer::brewer.pal(8, "Set1")
RColorBrewer::display.brewer.pal(4, "Set1")
cols <- cols[c(1:2, 4)]
names(cols) <- method_levels

position_dodge_width <- 0.2  # Adjust the width as needed

results_list <- list()
performance_list <- list()
sample_counts <- c(3, 5, 10, 20)

for (swap_no in 0:3)
{
  data_path <- paste0(validation_path, "/swap", swap_no, "/")
  results <- fread(paste0(data_path, "results.csv"))
  performance <- fread(paste0(data_path, "performance.csv"))

  results$celltype2 <- gsub(pattern = "_", replacement = " ", x = stringr::str_to_title(results$celltype))
  results$celltype2 <- gsub(pattern = "cd", replacement = "CD", results$celltype2)

  performance$celltype2 <- gsub(pattern = "_", replacement = " ", x = stringr::str_to_title(performance$celltype))
  performance$celltype2 <- gsub(pattern = "cd", replacement = "CD", performance$celltype2)

  results[,method := factor(method, levels=method_levels)]
  performance[,method := factor(method, levels=method_levels)]
  
  results[,n:=sample_counts[swap_no+1]]
  performance[,n:=sample_counts[swap_no+1]]
  
  results_list[[swap_no+1]] <- results
  performance_list[[swap_no+1]] <- performance
}

results_dt <- bind_rows(results_list)
performance_dt <- bind_rows(performance_list)

results_dt$n <- factor(results_dt$n, levels = sample_counts)
performance_dt$n <- factor(performance_dt$n, levels = sample_counts)

performance_dt$threshold <- factor(performance_dt$threshold, levels = c(0.01, 0.05, 0.1))

fwrite(results_dt, file = paste0(validation_path, "all_results.csv"))
fwrite(performance_dt, file = paste0(validation_path, "all_performances.csv"))

# Plot RMSE as n increases.
pl <- ggplot(results_dt, aes(n, rmse, fill=method)) + geom_boxplot() +
  scale_fill_manual(values = cols) + theme_bw() +
  theme(legend.title = element_blank())
pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "RMSE")
ggsave(filename = paste0(fig_path, "/RMSE.pdf"), pl)

# Plot CI coverage.
pl <- ggplot(results_dt, aes(n, coverage, fill=method)) + geom_boxplot() +
  scale_fill_manual(values = cols) + theme_bw() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 0.95, linetype = "dashed")
pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "Coverage")
ggsave(filename = paste0(fig_path, "/Coverage.pdf"), pl)

# Precision plot.
pl <- ggplot(performance_dt, aes(n, prec, fill=method)) + geom_boxplot() +
  scale_fill_manual(values = cols) + theme_bw() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = 0.95, linetype = "dashed")
pl <- pl + facet_grid(~ threshold)
pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "Precision")
ggsave(filename = paste0(fig_path, "/Precision.pdf"), pl, width = 18)

# False positive rates
pl <- ggplot(performance_dt, aes(n, fpr, fill=method)) + geom_boxplot() +
  scale_fill_manual(values = cols) + theme_bw() +
  theme(legend.title = element_blank())
pl <- pl + facet_grid(~ threshold)
pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "False positive rates")
ggsave(filename = paste0(fig_path, "/FPR.pdf"), pl, width = 18)

# False negative rates
pl <- ggplot(performance_dt, aes(n, fnr, fill=method)) + geom_boxplot() +
  scale_fill_manual(values = cols) + theme_cowplot() +
  theme(legend.title = element_blank())
pl <- pl + facet_grid(~ threshold)
pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "False negative rates")
ggsave(filename = paste0(fig_path, "/FNR.pdf"), pl, width = 18)

# Convey Precision, FNR, FPR for each method.
prec_dt <- performance_dt[,.(avg=mean(prec), std=sd(prec)),by=.(method, threshold, n)]
prec_dt[,ymin := avg - std]
prec_dt[,ymax := avg + std]
pl <- ggplot(prec_dt, aes(n, avg, col=threshold)) + 
  geom_line(aes(group=threshold), position = position_dodge(width=position_dodge_width)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width = 0.2, position = position_dodge(width=position_dodge_width)) +
  theme_cowplot() +
  #theme(legend.title = element_blank()) +
  labs(color="Threshold") +
  facet_grid(~method)
prec_pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "Precision")
ggsave(filename = paste0(fig_path, "/Precision2.pdf"), prec_pl, width = 18)

fnr_dt <- performance_dt[,.(avg=mean(fnr), std=sd(fnr)),by=.(method, threshold, n)]
fnr_dt[,ymin := avg - std]
fnr_dt[,ymax := avg + std]
pl <- ggplot(fnr_dt, aes(n, avg, col=threshold)) + 
  geom_line(aes(group=threshold), position = position_dodge(width=position_dodge_width)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width = 0.2, position = position_dodge(width=position_dodge_width)) +
  theme_cowplot() +
  #theme(legend.title = element_blank()) +
  labs(color="Threshold") +
  facet_grid(~method)
fnr_pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "FNR")
ggsave(filename = paste0(fig_path, "/FNR2.pdf"), fnr_pl, width = 18)

fpr_dt <- performance_dt[,.(avg=mean(fpr), std=sd(fpr)),by=.(method, threshold, n)]
fpr_dt[,ymin := avg - std]
fpr_dt[,ymax := avg + std]
pl <- ggplot(fpr_dt, aes(n, avg, col=threshold)) + 
  geom_line(aes(group=threshold), position = position_dodge(width=position_dodge_width)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width = 0.2, position = position_dodge(width=position_dodge_width)) +
  theme_cowplot() +
  #theme(legend.title = element_blank()) +
  labs(color="Threshold") +
  facet_grid(~method)
fpr_pl <- decorate_figure(pl, xlab_text = "Sample size", ylab_text = "FPR")
ggsave(filename = paste0(fig_path, "/FPR2.pdf"), fpr_pl, width = 18)

prec_pl <- prec_pl + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

fnr_pl <- fnr_pl + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank()
)

combined_plot <- plot_grid(prec_pl, fnr_pl, fpr_pl, ncol = 1, align = 'v', axis = 'tb')
gsave(filename = paste0(fig_path, "/Combined_Plot.pdf"), combined_plot, width = 18, height = 24) 

prec_dt[,metric:="Prec"]
fpr_dt[,metric:="FPR"]
fnr_dt[,metric:="FNR"]
metric_dt <- bind_rows(prec_dt, fpr_dt, fnr_dt)
pl <- ggplot(metric_dt, aes(n, avg, col=threshold)) +
  geom_point(position = position_dodge(width=position_dodge_width)) + 
  geom_line(aes(group=threshold), position = position_dodge(width=position_dodge_width)) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width = 0.2, position = position_dodge(width=position_dodge_width)) +
  theme_bw() +
  labs(x="Sample size", y = "Performance", color="Threshold") +
  facet_grid(metric~method, scales = "free_y") +
  theme(strip.text = element_text(size=20), legend.text = element_text(size=20)) + 
  theme(axis.title = element_text(size=20), legend.title = element_text(size=20))
ggsave(filename = paste0(fig_path, "/Combined_Plot.pdf"), pl, width = 18, height = 12)  

