library(tidyr)
library(ggplot2)
library(latex2exp)

save.results = function(results.mrse, results.err.phi=NULL, results.err.xi=NULL, save.path, filename) {
  results.mrse <- data.frame(do.call(rbind, lapply(results.mrse, function(r){do.call(rbind, r)})))
  df.mrse <- results.mrse %>% pivot_longer(1:length(methods), names_to = "method", values_to = "value")
  df.mrse$n <- as.factor(df.mrse$n)
  write.csv(df.mrse, paste0(save.path, "mrse.", filename, ".csv"), row.names = F)
  
  if (!is.null(results.err.phi)) {
    results.err.phi <- data.frame(do.call(rbind, lapply(results.err.phi, function(r){do.call(rbind, r)})))
    df.err.phi <- results.err.phi %>% pivot_longer(1:length(methods), names_to = "method", values_to = "value")
    df.err.phi$n <- as.factor(df.err.phi$n)
    write.csv(df.err.phi, paste0(save.path, "err.phi.", filename, ".csv"), row.names = F)
  }
  
  if (!is.null(results.err.xi)) {
    results.err.xi <- data.frame(do.call(rbind, lapply(results.err.xi, function(r){do.call(rbind, r)})))
    df.err.xi <- results.err.xi %>% pivot_longer(1:length(methods), names_to = "method", values_to = "value")
    df.err.xi$n <- as.factor(df.err.xi$n)
    write.csv(df.err.xi, paste0(save.path, "err.xi.", filename, ".csv"), row.names = F)
  }
}

generate.figures = function(tables.path, figures.path, filename, n.plot=NULL, only.mrse=F) {
  data <- read.csv(paste0(tables.path, "mrse.", filename, ".csv"))
  data$n <- as.factor(data$n)
  data$sparsity <- factor(data$sparsity, levels=c("Dense", "Low Sparsity", "Medium Sparsity", "High Sparsity"))

  ggplot(data, aes(x=n, y=value, fill=method)) + 
    geom_boxplot() + 
    xlab("N") + ylab("MRSE") + 
    scale_y_continuous(trans="log10") + 
    facet_grid(.~sparsity, scales="free") + 
    scale_fill_grey(start = 0.4, end=0.9, name="Method") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top") + 
    ggpubr::stat_compare_means(aes(group=method), label="p.signif", method="wilcox.test", paired=F, vjust = 2, hide.ns = T)
    
  
  ggsave(paste0(figures.path, "mrse.", filename, ".pdf"), width=11, height=3.5)
  
  if (!only.mrse) {
    data <- read.csv(paste0(tables.path, "err.phi.", filename, ".csv"))
    data$n <- as.factor(data$n)
    data$m <- as.factor(data$m)
    data$sparsity <- factor(data$sparsity, levels=c("Dense", "Low Sparsity", "Medium Sparsity", "High Sparsity"))
    
    ggplot(data[data$n == n.plot,], aes(x=m, y=value, fill=method)) + 
      geom_boxplot() + 
      xlab("m") + ylab(TeX(r"($MSE(\hat{\phi}_i)$)")) + 
      scale_y_continuous(trans="log10") + 
      facet_grid(.~sparsity) + 
      scale_fill_grey(start = 0.4, end=0.9) + 
      theme_bw(base_size = 14) +
      theme(legend.position = "top")
    
    ggsave(paste0(figures.path, "err.phi.n.", n.plot, ".", filename, ".pdf"), width=11, height=3.5)
    
    data <- read.csv(paste0(tables.path, "err.xi.", filename, ".csv"))
    data$n <- as.factor(data$n)
    data$m <- as.factor(data$m)
    data$sparsity <- factor(data$sparsity, levels=c("Dense", "Low Sparsity", "Medium Sparsity", "High Sparsity"))
    
    ggplot(data[data$n == n.plot,], aes(x=m, y=value, fill=method)) + 
      geom_boxplot() + 
      xlab("m") + ylab(TeX(r"($MSE(\hat{\xi}_i)$)")) + 
      scale_y_continuous(trans="log10") + 
      facet_grid(.~sparsity) + 
      scale_fill_grey(start = 0.4, end=0.9) + 
      theme_bw(base_size = 14) +
      theme(legend.position = "top")
    
    ggsave(paste0(figures.path, "err.xi.n.", n.plot, ".", filename, ".pdf"), width=11, height=3.5)
  }
}
