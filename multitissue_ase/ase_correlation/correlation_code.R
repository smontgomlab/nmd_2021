library(magrittr)
library(scales)
library(openxlsx)
library(tidyverse)

tissue_cors_df <- readRDS("./tissue_cors_df_new.rda")

tissue_cors_df %<>% separate(Comb, into = c("Tissue1", "sep", "Tissue2"), sep = "_") %>% select(-sep)
tissue_cors_df$cor %<>% signif(digits = 3)

tissue_cors_self <- tibble(cor = NA, p = NA, Z = NA, t = NA, nObs = NA, Variant = NA, Tissue1 = unique(c(tissue_cors_df$Tissue1, tissue_cors_df$Tissue2)), Tissue2 = unique(c(tissue_cors_df$Tissue1, tissue_cors_df$Tissue2)))

tissue_cors_stop_gained <- filter(tissue_cors_df, Variant == "stop_gained")
tissue_cors_missense <- filter(tissue_cors_df, Variant == "missense_variant")
colnames(tissue_cors_missense)[c(7,8)] <- c("Tissue2", "Tissue1")
tissue_cors_missense %<>% select(cor:Variant, Tissue1, Tissue2)

tissue_cors_stop_vs_missense <- bind_rows(tissue_cors_stop_gained, tissue_cors_missense, tissue_cors_self)
tissue_cors_sm_plot <- ggplot(tissue_cors_stop_vs_missense, aes(Tissue1, Tissue2, fill = cor, label = cor)) +
    geom_tile() +
    geom_text() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_fill_continuous(low = "white", high = muted("green"), na.value = "black", limits = c(0,1))

CairoPDF("cor_stop_vs_missense", height = 30, width = 30)
    print(tissue_cors_sm_plot)
dev.off()

tissue_cors_synonymous <- filter(tissue_cors_df, Variant == "synonymous_variant")
colnames(tissue_cors_synonymous)[c(7,8)] <- c("Tissue2", "Tissue1")
tissue_cors_synonymous %<>% select(cor:Variant, Tissue1, Tissue2)

tissue_cors_stop_vs_synonymous <- bind_rows(tissue_cors_stop_gained, tissue_cors_synonymous, tissue_cors_self)
tissue_cors_sm_plot <- ggplot(tissue_cors_stop_vs_synonymous, aes(Tissue1, Tissue2, fill = cor, label = cor)) +
    geom_tile() +
    geom_text() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_fill_continuous(low = "white", high = muted("green"), na.value = "black", limits = c(0,1))

CairoPDF("cor_stop_vs_synonymous", height = 30, width = 30)
    print(tissue_cors_sm_plot)
dev.off()

tissue_cors_intron <- filter(tissue_cors_df, Variant == "intron_variant")
colnames(tissue_cors_intron)[c(7,8)] <- c("Tissue2", "Tissue1")
tissue_cors_intron %<>% select(cor:Variant, Tissue1, Tissue2)

tissue_cors_stop_vs_intron <- bind_rows(tissue_cors_stop_gained, tissue_cors_intron, tissue_cors_self)
tissue_cors_sm_plot <- ggplot(tissue_cors_stop_vs_intron, aes(Tissue1, Tissue2, fill = cor, label = cor)) +
    geom_tile() +
    geom_text() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_fill_continuous(low = "white", high = muted("green"), na.value = "black", limits = c(0,1))

CairoPDF("cor_stop_vs_intron", height = 30, width = 30)
    print(tissue_cors_sm_plot)
dev.off()

tissue_cors_3_prime_UTR <- filter(tissue_cors_df, Variant == "3_prime_UTR_variant")
colnames(tissue_cors_3_prime_UTR)[c(7,8)] <- c("Tissue2", "Tissue1")
tissue_cors_3_prime_UTR %<>% select(cor:Variant, Tissue1, Tissue2)

tissue_cors_stop_vs_3_prime_UTR <- bind_rows(tissue_cors_stop_gained, tissue_cors_3_prime_UTR, tissue_cors_self)
tissue_cors_sm_plot <- ggplot(tissue_cors_stop_vs_3_prime_UTR, aes(Tissue1, Tissue2, fill = cor, label = cor)) +
    geom_tile() +
    geom_text() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_fill_continuous(low = "white", high = muted("green"), na.value = "black", limits = c(0,1))

CairoPDF("cor_stop_vs_3_prime_UTR", height = 30, width = 30)
    print(tissue_cors_sm_plot)
dev.off()

tissue_cors_non_coding_transcript_exon <- filter(tissue_cors_df, Variant == "non_coding_transcript_exon_variant")
colnames(tissue_cors_non_coding_transcript_exon)[c(7,8)] <- c("Tissue2", "Tissue1")
tissue_cors_non_coding_transcript_exon %<>% select(cor:Variant, Tissue1, Tissue2)

tissue_cors_stop_vs_non_coding_transcript_exon <- bind_rows(tissue_cors_stop_gained, tissue_cors_non_coding_transcript_exon, tissue_cors_self)
tissue_cors_sm_plot <- ggplot(tissue_cors_stop_vs_non_coding_transcript_exon, aes(Tissue1, Tissue2, fill = cor, label = cor)) +
    geom_tile() +
    geom_text() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_fill_continuous(low = "white", high = muted("green"), na.value = "black", limits = c(0,1))

CairoPDF("cor_stop_vs_non_coding_transcript_exon", height = 30, width = 30)
    print(tissue_cors_sm_plot)
dev.off()

AddWorkbook <- function(index, tissue_cors_list, wb) {
    sheet_name <- names(tissue_cors_list)[index]
    tissue_cor_table <- tissue_cors_list[[index]]
    tissue_cor_table %<>% mutate(across(where(is.numeric)))

    addWorksheet(wb = wb, sheetName = sheet_name)
    writeDataTable(wb = wb, sheet = index, x = tissue_cor_table) 
    setColWidths(wb, index, cols = 1:ncol(tissue_cor_table), widths = "auto") 
}

tissue_cors_split <- split(tissue_cors_df, tissue_cors_df$Variant) 
names(tissue_cors_split) %<>% str_replace_all("non_coding", "nc")

wb <- createWorkbook()
map(1:length(tissue_cors_split), AddWorkbook, tissue_cors_split, wb)
saveWorkbook(wb, "tissue_cors_supp.xlsx", overwrite = TRUE)

