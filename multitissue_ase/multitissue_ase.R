library(parallel)
library(pbmcapply)
library(ComplexHeatmap)

library(ggdendro)
library(Cairo)
library(magrittr)
library(tidyverse)

AggregateTissues <- function(wasp_rows) {
    tissue_id <- unique(wasp_rows$TISSUE_ID)
    count_cols <- select(wasp_rows, REF_COUNT:TOTAL_COUNT) %>% colSums %>% as.matrix %>% t %>% data.frame
    ref_ratio <- count_cols$REF_COUNT / count_cols$TOTAL_COUNT
    annot_info <- select(wasp_rows, CHR:ALT_ALLELE) %>% distinct
    annot_info$TISSUE_ID <- tissue_id
    annot_tibble <- cbind(annot_info, count_cols) %>% as_tibble
    annot_tibble$REF_RATIO <- ref_ratio
    annot_tibble$Count <- nrow(wasp_rows)
    annot_tibble
}

AggregateVariant <- function(variant_id, wasp_join) {
    wasp_variant <- filter(wasp_join, VARIANT_ID == variant_id)
    variant_aggregated <- split(wasp_variant, wasp_variant$TISSUE_ID) %>% map_df(AggregateTissues)
    variant_aggregated
}

variants_annot <- readRDS("../annotate_variants/variants_annotated.rda")
variants_annot$VARIANT_ID <- str_c("chr", variants_annot$Chrom, "_", variants_annot$Pos, "_", variants_annot$Ref, "_", variants_annot$Alt, "_b38")
variants_dup <- group_by(variants_annot, VARIANT_ID, Gene) %>% 
    summarise(Count = n()) %>% 
    group_by(VARIANT_ID) %>% 
    summarise(Count = n()) %>% 
    filter(Count > 1)

variants_rare <- filter(variants_annot, gnomAD_AF < 0.01 & !is_in(VARIANT_ID, variants_dup$VARIANT_ID))

wasp_in <- read_tsv("./all_subjects_WASP_stopgain.tsv")

sample_attributes <- read_tsv("./GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt") %>%
    filter(ANALYTE_TYPE == "RNA:Total RNA") %>%
    select(SAMPID, SMRIN, SMTSISCH, SMGEBTCH)
colnames(sample_attributes)[1] <- "SAMPLE_ID"
subject_attributes <- read_tsv("./GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt") %>%
    select(SUBJID, SEX, AGE, DTHHRDY)
colnames(subject_attributes)[1] <- "SUBJECT_ID"

#genotypes <- readRDS("./main_data_genotype.rda") %>% select(SAMPLE_ID, TISSUE_ID, Genotype)
#genotype_split <- str_split_fixed(genotypes$Genotype, ",", 4)[,1] 
#genotypes$Genotype_min <- str_split_fixed(genotype_split, ":", 2)[,1]
#genotypes$Het_Status <- genotypes$Genotype_min == "0/1"

wasp_annot <- left_join(wasp_in, sample_attributes) %>% left_join(subject_attributes)

wasp_rare <- filter(wasp_annot, is_in(VARIANT_ID, variants_rare$VARIANT_ID)) 
wasp_rare_filter <- filter(wasp_rare, !is_in(TISSUE_ID, c("KDNMDL", "FLLPNT", "CVSEND", "CVXECT", "BLDDER")))
wasp_rare_filter$split_var <- str_c(wasp_rare_filter$VARIANT_ID, wasp_rare_filter$SUBJECT_ID, sep = "_")

wasp_split <- split(wasp_rare_filter, wasp_rare_filter$split_var)
wasp_split_nrows <- map_int(wasp_split, nrow)
wasp_split_multi <- wasp_split[wasp_split_nrows > 1]
wasp_split_counts <- map(wasp_split, select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)

#Unaggregated
source("ASE_29Dec2014_new.R")
all_ase <- pbmclapply(wasp_split_counts, gtm, model.strong.ase = FALSE, mc.cores = 8)
saveRDS(all_ase, "./all_ase.rda")

all_ase <- readRDS("../../NMD_old/NMD4/all_ase.rda")
ase_posteriors <- map(all_ase, extract2, "indiv.posteriors") %>% 
    map(magrittr::extract, "MODASE", TRUE) %>% 
    map(t) %>% reshape2::melt() %>% select(-Var2)
colnames(ase_posteriors) <- c("TISSUE_ID", "MODASE", "split_var")

wasp_modase_join <- inner_join(wasp_rare_filter, ase_posteriors) 
wasp_modase_join$MODASE_cat <- as.integer(wasp_modase_join$MODASE >= 0.8 & wasp_modase_join$REF_RATIO > 0.5) %>% factor %>% fct_recode(imbalanced = "1", balanced = "0")
wasp_modase_join$imbalanced_cat <- as.integer(wasp_modase_join$REF_RATIO >= 0.65) %>% factor %>% fct_recode(imbalanced = "1", balanced = "0")

wasp_final_noagg <- left_join(wasp_modase_join, variants_rare)
write_tsv(wasp_final_noagg, "ase_unaggegated.txt")

#Aggregated
wasp_annot2 <- readRDS("../predictive_models/training_data.rda")
wasp_filter <- filter(wasp_annot2, !is_in(TISSUE_ID, c("KDNMDL", "FLLPNT", "CVSEND", "CVXECT", "BLDDER")))
wasp_male <- filter(wasp_filter, SEX == 1 & !str_detect(TISSUE_ID, "BRN"))
wasp_female <- filter(wasp_filter, SEX == 2 & !str_detect(TISSUE_ID, "BRN"))
wasp_brain <- filter(wasp_filter, str_detect(TISSUE_ID, "BRN"))

wasp_group_brain <- group_by(wasp_brain, SUBJECT_ID, VARIANT_ID) %>% summarise(Count = n()) %>% group_by(VARIANT_ID) %>% summarise(Count = n()) 

wasp_multi_group_brain <- filter(wasp_group_brain, Count > 1) %>% arrange(desc(Count))
wasp_multi_agg_brain <- map_df(wasp_multi_group_brain$VARIANT_ID, AggregateVariant, wasp_brain)

wasp_multi_counts_brain <- split(wasp_multi_agg_brain, wasp_multi_agg_brain$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_multi_lengths_brain <- map_int(wasp_multi_counts_brain, nrow)
wasp_multi_multi_brain <- wasp_multi_counts_brain[wasp_multi_lengths_brain > 1]

#multi_ase_brain <- pbmclapply(wasp_multi_multi_brain, gtm, model.strong.ase = FALSE, mc.cores = 8)
#saveRDS(multi_ase, "./multi_ase.rda")

#Extract singleton variants from unaggregated data
wasp_unique_group_brain <- filter(wasp_group_brain, Count == 1) %>% arrange(desc(Count))

wasp_unique_filter_brain <- filter(wasp_brain, is_in(VARIANT_ID, wasp_unique_group_brain$VARIANT_ID))
wasp_unique_counts_brain <- split(wasp_unique_filter_brain, wasp_unique_filter_brain$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_unique_lengths_brain <- map_int(wasp_unique_counts_brain, nrow)
wasp_unique_multi_brain <- wasp_unique_counts_brain[wasp_unique_lengths_brain > 1]

#Combine single variants and aggregated variants in multiple subjects
wasp_agg_counts_brain <- c(wasp_unique_multi_brain, wasp_multi_multi_brain)

gtm_star_output_brain <- gtm.star(wasp_agg_counts_brain, model.strong.ase = FALSE)
saveRDS(gtm_star_output_brain, "./gtm_star_output_brain.rda")

brain_clust_plot <- as.dist(gtm_star_output_brain$distances) %>% hclust(method = "average") 

CairoPDF("brain_clust", bg = "transparent")
     plot(brain_clust_plot, sub = "", xlab = "", ylab = "distance", main = "")
dev.off()

brain_mds <- as.dist(gtm_star_output_brain$distances) %>% cmdscale %>%
    as.data.frame %>% rownames_to_column("Tissue")
brain_mds_plot <- ggplot(brain_mds, aes(V1, V2, label = Tissue)) +
    geom_text() +
    theme_classic() +
    xlab("PC1") +
    ylab("PC2") +
    theme(panel.border = element_rect(fill = NA),
        plot.background = element_blank())

CairoPDF("brain_mds", bg = "transparent")
    print(brain_mds_plot)
dev.off()

#Male
wasp_group_male <- group_by(wasp_male, SUBJECT_ID, VARIANT_ID) %>% summarise(Count = n()) %>% group_by(VARIANT_ID) %>% summarise(Count = n()) 

wasp_multi_group_male <- filter(wasp_group_male, Count > 1) %>% arrange(desc(Count))
wasp_multi_agg_male <- map_df(wasp_multi_group_male$VARIANT_ID, AggregateVariant, wasp_male)

wasp_multi_counts_male <- split(wasp_multi_agg_male, wasp_multi_agg_male$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_multi_lengths_male <- map_int(wasp_multi_counts_male, nrow)
wasp_multi_multi_male <- wasp_multi_counts_male[wasp_multi_lengths_male > 1]

#multi_ase_male <- pbmclapply(wasp_multi_multi_male, gtm, model.strong.ase = FALSE, mc.cores = 8)
#saveRDS(multi_ase, "./multi_ase.rda")

#Extract singleton variants from unaggregated data
wasp_unique_group_male <- filter(wasp_group_male, Count == 1) %>% arrange(desc(Count))

wasp_unique_filter_male <- filter(wasp_male, is_in(VARIANT_ID, wasp_unique_group_male$VARIANT_ID))
wasp_unique_counts_male <- split(wasp_unique_filter_male, wasp_unique_filter_male$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_unique_lengths_male <- map_int(wasp_unique_counts_male, nrow)
wasp_unique_multi_male <- wasp_unique_counts_male[wasp_unique_lengths_male > 1]

#Combine single variants and aggregated variants in multiple subjects
wasp_agg_counts_male <- c(wasp_unique_multi_male, wasp_multi_multi_male)

gtm_star_output_male <- gtm.star(wasp_agg_counts_male, model.strong.ase = FALSE)
saveRDS(gtm_star_output_male, "./gtm_star_output_male.rda")

male_clust_plot <- as.dist(gtm_star_output_male$distances) %>% hclust(method = "average") 

CairoPDF("male_clust", bg = "transparent")
     plot(male_clust_plot, sub = "", xlab = "", ylab = "distance", main = "")
dev.off()

male_mds <- as.dist(gtm_star_output_male$distances) %>% cmdscale %>%
    as.data.frame %>% rownames_to_column("Tissue")
male_mds_plot <- ggplot(male_mds, aes(V1, V2, label = Tissue)) +
    geom_text() +
    theme_classic() +
    xlab("PC1") +
    ylab("PC2") +
    theme(panel.border = element_rect(fill = NA),
        plot.background = element_blank())

CairoPDF("male_mds", bg = "transparent")
    print(male_mds_plot)
dev.off()

#Female
wasp_group_female <- group_by(wasp_female, SUBJECT_ID, VARIANT_ID) %>% summarise(Count = n()) %>% group_by(VARIANT_ID) %>% summarise(Count = n()) 

wasp_multi_group_female <- filter(wasp_group_female, Count > 1) %>% arrange(desc(Count))
wasp_multi_agg_female <- map_df(wasp_multi_group_female$VARIANT_ID, AggregateVariant, wasp_female)

wasp_multi_counts_female <- split(wasp_multi_agg_female, wasp_multi_agg_female$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_multi_lengths_female <- map_int(wasp_multi_counts_female, nrow)
wasp_multi_multi_female <- wasp_multi_counts_female[wasp_multi_lengths_female > 1]

#multi_ase_female <- pbmclapply(wasp_multi_multi_female, gtm, model.strong.ase = FALSE, mc.cores = 8)
#saveRDS(multi_ase, "./multi_ase.rda")

#Extract singleton variants from unaggregated data
wasp_unique_group_female <- filter(wasp_group_female, Count == 1) %>% arrange(desc(Count))

wasp_unique_filter_female <- filter(wasp_female, is_in(VARIANT_ID, wasp_unique_group_female$VARIANT_ID))
wasp_unique_counts_female <- split(wasp_unique_filter_female, wasp_unique_filter_female$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_unique_lengths_female <- map_int(wasp_unique_counts_female, nrow)
wasp_unique_multi_female <- wasp_unique_counts_female[wasp_unique_lengths_female > 1]

#Combine single variants and aggregated variants in multiple subjects
wasp_agg_counts_female <- c(wasp_unique_multi_female, wasp_multi_multi_female)

gtm_star_output_female <- gtm.star(wasp_agg_counts_female, model.strong.ase = FALSE)
saveRDS(gtm_star_output_female, "./gtm_star_output_female.rda")
gtm_star_output_female <- readRDS("./gtm_star_output_female.rda")

female_clust_plot <- as.dist(gtm_star_output_female$distances) %>% hclust(method = "average") 

CairoPDF("female_clust", bg = "transparent")
     plot(female_clust_plot, sub = "", xlab = "", ylab = "distance", main = "")
dev.off()

female_mds <- as.dist(gtm_star_output_female$distances) %>% cmdscale %>%
    as.data.frame %>% rownames_to_column("Tissue")
female_mds_plot <- ggplot(female_mds, aes(V1, V2, label = Tissue)) +
    geom_text() +
    theme_classic() +
    xlab("PC1") +
    ylab("PC2") +
    theme(panel.border = element_rect(fill = NA),
        plot.background = element_blank())

CairoPDF("female_mds", bg = "transparent")
    print(female_mds_plot)
dev.off()

#No sex organs or brain
wasp_filter_strict <- filter(wasp_filter, !is_in(TISSUE_ID, c("TESTIS", "PRSTTE", "VAGINA", "UTERUS", "OVARY", "PTTARY")) & 
    !str_detect(TISSUE_ID, "BRN"))

wasp_group_strict <- group_by(wasp_filter_strict, SUBJECT_ID, VARIANT_ID) %>% summarise(Count = n()) %>% group_by(VARIANT_ID) %>% summarise(Count = n()) 
wasp_subject_strict <- group_by(wasp_filter_strict, SUBJECT_ID, TISSUE_ID) %>% tally %>% group_by(SUBJECT_ID) %>% tally

subject_plot <- ggplot(wasp_subject_strict, aes(n)) +
    geom_histogram(stat = "bin", binwidth = 1, color = "black", fill = "white") +
    theme_classic() +
    scale_x_continuous(breaks = seq(1,28, by = 1))

CairoPDF("subject_hist", width = 10)
    print(subject_plot)
dev.off()

wasp_multi_group_strict <- filter(wasp_group_strict, Count > 1) %>% arrange(desc(Count))
wasp_multi_agg_strict <- map_df(wasp_multi_group_strict$VARIANT_ID, AggregateVariant, wasp_filter_strict)

wasp_multi_counts_strict <- split(wasp_multi_agg_strict, wasp_multi_agg_strict$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_multi_lengths_strict <- map_int(wasp_multi_counts_strict, nrow)
wasp_multi_multi_strict <- wasp_multi_counts_strict[wasp_multi_lengths_strict > 1]

#Extract singleton variants from unaggregated data
wasp_unique_group_strict <- filter(wasp_group_strict, Count == 1) %>% arrange(desc(Count))

wasp_unique_filter_strict <- filter(wasp_filter_strict, is_in(VARIANT_ID, wasp_unique_group_strict$VARIANT_ID))
wasp_unique_counts_strict <- split(wasp_unique_filter_strict, wasp_unique_filter_strict$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_unique_lengths_strict <- map_int(wasp_unique_counts_strict, nrow)
wasp_unique_multi_strict <- wasp_unique_counts_strict[wasp_unique_lengths_strict > 1]

#Combine single variants and aggregated variants in multiple subjects
wasp_agg_counts_strict <- c(wasp_unique_multi_strict, wasp_multi_multi_strict)

gtm_star_output_strict <- gtm.star(wasp_agg_counts_strict, model.strong.ase = FALSE)
saveRDS(gtm_star_output_strict, "./gtm_star_output_nosex_nobrain.rda")
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

#Random split of data
subject_id <- unique(wasp_filter_strict$SUBJECT_ID)
set.seed(12345)
subject_id_half <- sample(subject_id, size = ceiling(length(subject_id) / 2))

wasp_filter_strict1 <- filter(wasp_filter_strict, is_in(SUBJECT_ID, subject_id_half))
wasp_filter_strict2 <- filter(wasp_filter_strict, !is_in(SUBJECT_ID, subject_id_half))

SharingHeatmap <- function(data_subset, filename, height = 8, width = 6) {
    subject_sharing_matrix <- select(data_subset, SUBJECT_ID, TISSUE_ID) %>%
        distinct %>%
        mutate(Count = 1) %>%
        pivot_wider(names_from = TISSUE_ID, values_from = Count, values_fill = list(Count = 0)) %>%
        select(-SUBJECT_ID) %>%
        mutate_all(as.integer) %>%
        as.matrix

    print(dim(subject_sharing_matrix))
    print(typeof(subject_sharing_matrix))

    plot_heatmap <- Heatmap(subject_sharing_matrix, 
                show_row_dend = FALSE, 
                col = c("white", "black"),
                show_heatmap_legend = FALSE)
    plot_heatmap

    #CairoPDF(filename, height = height, width = width, bg = "transparent")
        #print(plot_heatmap)
    #dev.off()
}

SharingHeatmap2 <- function(data_subset, filename, height = 8, width = 6) {
    subject_sharing_matrix <- select(data_subset, `#Uploaded_variation`, TISSUE_ID) %>%
        distinct %>%
        mutate(Count = 1) %>%
        pivot_wider(names_from = TISSUE_ID, values_from = Count, values_fill = list(Count = 0)) %>%
        select(-`#Uploaded_variation`) %>%
        mutate_all(as.integer) %>%
        as.matrix

    plot_heatmap <- Heatmap(subject_sharing_matrix, 
        use_raster = F,
        show_row_dend = FALSE, 
        col = c("white", "black"),
        show_heatmap_legend = FALSE)

    CairoPDF(filename, height = height, width = width, bg = "transparent")
        print(plot_heatmap)
    dev.off()
}

SharingHeatmap(wasp_filter_strict, "subject_sharing_nosex_nobrain")
SharingHeatmap(wasp_filter_strict1, "subject_sharing_nosex_nobrain_random_subset1")
SharingHeatmap(wasp_filter_strict2, "subject_sharing_nosex_nobrain_random_subset2")

wasp_filter_strict_male <- filter(wasp_filter_strict, SEX == 1)
wasp_filter_strict_female <- filter(wasp_filter_strict, SEX == 2)

SharingHeatmap(wasp_filter_strict_male, "subject_sharing_nosex_nobrain_male")
SharingHeatmap(wasp_filter_strict_female, "subject_sharing_nosex_nobrain_female")

SharingHeatmap(wasp_filter, "subject_sharing_full", height = 8, width = 8)
wasp_male2 <- filter(wasp_filter, SEX == 1)
wasp_female2 <- filter(wasp_filter, SEX == 2)
SharingHeatmap(wasp_male2, "subject_sharing_male", width = 7)
SharingHeatmap(wasp_female2, "subject_sharing_female", width = 7)

wasp_nobrain <- filter(wasp_filter, !str_detect(TISSUE_ID, "BRN"))
SharingHeatmap(wasp_nobrain, "subject_sharing_nobrain")
SharingHeatmap(wasp_male, "subject_sharing_nobrain_male")
SharingHeatmap(wasp_female, "subject_sharing_nobrain_female")

wasp_brain <- filter(wasp_filter, str_detect(TISSUE_ID, "BRN"))
SharingHeatmap(wasp_brain, "subject_sharing_brain_only")

SharingHeatmap2(wasp_filter, "variant_sharing_all", height = 12, width = 10)
test2 <- SharingHeatmap(wasp_filter, "")

print(test1)
print(test2)

plot_heatmap <- Heatmap(sharing_matrix, 
                        use_raster = F,
    show_row_dend = FALSE, 
    col = c("white", "black"),
    show_heatmap_legend = FALSE)


STOP
#subject_sharing_clust <- dist(t(subject_sharing_matrix)) %>% hclust

#CairoPDF("subject_sharing_hclust")
    #plot(subject_sharing_clust)
#dev.off()

#subject_sharing_plot <- ggplot(wasp_filter_strict, aes(TISSUE_ID, SUBJECT_ID)) +
    #geom_tile() +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_blank()
    #)

#subject_sharing_plot1 <- ggplot(wasp_filter_strict1, aes(TISSUE_ID, SUBJECT_ID)) +
    #geom_tile() +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_blank()
    #)

#CairoPDF("subject_sharing1", height = 8, width = 6)
    #print(subject_sharing_plot1)
#dev.off()

#subject_sharing_plot2 <- ggplot(wasp_filter_strict2, aes(TISSUE_ID, SUBJECT_ID)) +
    #geom_tile() +
    #theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_blank()
    #)

#CairoPDF("subject_sharing2", height = 8, width = 6)
    #print(subject_sharing_plot2)
#dev.off()

wasp_group_strict1 <- group_by(wasp_filter_strict1, SUBJECT_ID, VARIANT_ID) %>% summarise(Count = n()) %>% group_by(VARIANT_ID) %>% summarise(Count = n()) 

wasp_multi_group_strict1 <- filter(wasp_group_strict1, Count > 1) %>% arrange(desc(Count))
wasp_multi_agg_strict1 <- map_df(wasp_multi_group_strict1$VARIANT_ID, AggregateVariant, wasp_filter_strict1)

wasp_multi_counts_strict1 <- split(wasp_multi_agg_strict1, wasp_multi_agg_strict1$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_multi_lengths_strict1 <- map_int(wasp_multi_counts_strict1, nrow)
wasp_multi_multi_strict1 <- wasp_multi_counts_strict1[wasp_multi_lengths_strict1 > 1]

#Extract singleton variants from unaggregated data
wasp_unique_group_strict1 <- filter(wasp_group_strict1, Count == 1) %>% arrange(desc(Count))

wasp_unique_filter_strict1 <- filter(wasp_filter_strict1, is_in(VARIANT_ID, wasp_unique_group_strict1$VARIANT_ID))
wasp_unique_counts_strict1 <- split(wasp_unique_filter_strict1, wasp_unique_filter_strict1$VARIANT_ID) %>% map(select, TISSUE_ID, REF_COUNT, ALT_COUNT) %>% map(column_to_rownames, "TISSUE_ID") %>% map(as.matrix)
wasp_unique_lengths_strict1 <- map_int(wasp_unique_counts_strict1, nrow)
wasp_unique_multi_strict1 <- wasp_unique_counts_strict1[wasp_unique_lengths_strict1 > 1]

#Combine single variants and aggregated variants in multiple subjects
wasp_agg_counts_strict1 <- c(wasp_unique_multi_strict1, wasp_multi_multi_strict1)

gtm_star_output_strict1 <- gtm.star(wasp_agg_counts_strict1, model.strong.ase = FALSE)
saveRDS(gtm_star_output_strict1, "./gtm_star_output_strict1.rda")
#gtm_star_output_strict1 <- readRDS("./gtm_star_output_strict1.rda")

strict1_clust_plot <- as.dist(gtm_star_output_strict1$distances) %>% hclust(method = "average") 

CairoPDF("strict1_clust", bg = "transparent")
     plot(strict1_clust_plot, sub = "", xlab = "", ylab = "distance", main = "")
dev.off()

strict1_mds <- as.dist(gtm_star_output_strict1$distances) %>% cmdscale %>%
    as.data.frame %>% rownames_to_column("Tissue")
strict1_mds_plot <- ggplot(strict1_mds, aes(V1, V2, label = Tissue)) +
    geom_text() +
    theme_classic() +
    xlab("PC1") +
    ylab("PC2") +
    theme(panel.border = element_rect(fill = NA),
        plot.background = element_blank())

CairoPDF("strict1_mds", bg = "transparent")
    print(strict1_mds_plot)
dev.off()
