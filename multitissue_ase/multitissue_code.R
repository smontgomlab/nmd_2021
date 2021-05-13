library(ggdendro)
library(Cairo)
library(magrittr)
library(tidyverse)

gtm_star_output_nosex_nobrain <- readRDS("./gtm_star_output_nosex_nobrain.rda")

nosex_nobrain_clust <- as.dist(gtm_star_output_nosex_nobrain$distances) %>% hclust(method = "average") 
nosex_nobrain_dendro <- dendro_data(nosex_nobrain_clust)
nosex_nobrain_clust_plot <- ggplot() +
    geom_segment(data = nosex_nobrain_dendro$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(breaks = seq_along(nosex_nobrain_dendro$labels$label), labels = nosex_nobrain_dendro$labels$label) +
    ylab("Height") +
    theme_classic() +
    theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, margin = margin(t = -15)),
        plot.title = element_blank(),
        axis.title.x = element_blank()) 

CairoPDF("nosex_nobrain_clust", bg = "transparent", height = 6, width = 10)
     print(nosex_nobrain_clust_plot)
dev.off()

nosex_nobrain_mds <- as.dist(gtm_star_output_nosex_nobrain$distances) %>% cmdscale %>%
    as.data.frame %>% rownames_to_column("Tissue")
nosex_nobrain_mds_plot <- ggplot(nosex_nobrain_mds, aes(V1, V2, label = Tissue)) +
    geom_text() +
    theme_classic() +
    xlab("PC1") +
    ylab("PC2") +
    theme(panel.border = element_rect(fill = NA),
        plot.background = element_blank())

CairoPDF("nosex_nobrain_mds", bg = "transparent")
    print(nosex_nobrain_mds_plot)
dev.off()

