# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(geomtextpath))

# Load data
kmers.df <- read_tsv("hist.txt", col_names = c("depth","count","log_scale", "gc"), comment = '#', show_col_types = FALSE)
peaks.df <- read_tsv("peaks.txt", col_names = c("start","center","stop", "max", "volume", "gc"), comment = '#', show_col_types = FALSE)
peaks1 <- as.numeric(peaks.df[1,2])
peaks2 <- as.numeric(peaks.df[2,2])

# Make plot
ggplot(data=kmers.df, aes(x=depth, y=log_scale)) +
    #geom_vline(data=peaks.df, mapping=aes(xintercept=center), color="grey") +
    geom_textvline(label = paste("Peak 1:", peaks1), xintercept = peaks1, color = "orchid", linewidth = 1) +
    geom_textvline(label = paste("Peak 2:", peaks2), xintercept = peaks2, color = "orchid", linewidth = 1) +
    geom_path() +
    scale_x_continuous(limits = c(0,500), expand = expansion(add=5)) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = "23-mer frequency", y = "count (log scale)") +
    theme_classic() + theme()
ggsave("hist.pdf", height=6, width=8)

# Add GC content tracer
kmers.df$gc_mod <- (kmers.df$gc/100) * (max(kmers.df[c(0:500),]$log_scale))
ggplot(data=kmers.df, aes(x=depth, y=log_scale)) +
    geom_hline(yintercept=max(kmers.df[c(0:500),]$log_scale)*0.5, linetype="dashed", color="orchid", alpha = 0.2) +
    geom_step(aes(x=depth, y=gc_mod), color = "orchid", alpha = 0.3) +
    #geom_vline(data=peaks.df, mapping=aes(xintercept=center), color="grey") +
    geom_textvline(label = paste("Peak 1:", peaks1), xintercept = peaks1, color = "orchid", linewidth = 1) +
    geom_textvline(label = paste("Peak 2:", peaks2), xintercept = peaks2, color = "orchid", linewidth = 1) +
    geom_path() +
    scale_x_continuous(limits = c(0,500), expand = expansion(add=5)) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x = "23-mer frequency", y = "count (log scale)") +
    theme_classic() + theme()
ggsave("hist_gc.pdf", height=6, width=8)
