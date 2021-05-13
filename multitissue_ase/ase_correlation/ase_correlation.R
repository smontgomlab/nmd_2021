library(WGCNA)
#library(parallel)
library(magrittr)
library(scales)
library(rtracklayer)
library(tidyverse)

all_annot <- read_tsv("../../annotate_variants/variants_unambiguous.txt", col_types = cols(.default = col_character()))
colnames(all_annot)[1] <- "Variant_ID"
all_annot_fix <- separate(all_annot, Variant_ID, into = c("CHR", "POS", "GT"), sep = "_", remove = FALSE) %>%
    separate(GT, c("REF_ALLELE", "ALT_ALLELE"))
all_annot_fix$CHR <- str_c("chr", all_annot_fix$CHR)
all_annot_fix$POS %<>% as.integer
all_annot_reduce <- select(all_annot_fix, CHR:ALT_ALLELE, Variant_ID, Consequence)
saveRDS(all_annot_reduce, "./all_annot_reduce.rda")
all_annot_reduce <- readRDS("all_annot_reduce.rda")

gnomad_annot <- read_tsv("../../annotate_variants/variants_gnomad_only.txt", col_types = cols(.default = col_character())) %>% distinct
gnomad_annot$gnomAD_AF %<>% str_replace_all("\\-", "0") %>% as.numeric
gnomad_annot_maf01 <- filter(gnomad_annot, gnomAD_AF < 0.01) %>% distinct

STOP
#test_read$POS %<>% as.character

#test_join <- left_join(test_read, all_annot_fix)

#test_agree <- filter(test_join, VARIANT_ANNOTATION == Consequence)
#test_stopgain <- filter(test_join, VARIANT_ANNOTATION == "stop_gained" | Consequence == "stop_gained")

#discordant <- filter(test_join, VARIANT_ANNOTATION != Consequence)

#SplitTissues <- function(tissue_rows, tissue_id, subject_id) {
    #write_tsv(tissue_rows, str_c("../raw_ase", subject_id, tissue_id, sep = "/"))
#}

#MapSplitTissues <- function(sample_file) {
    #sample_read <- read_tsv(sample_file)
    #subject_id <- unique(sample_read$SUBJECT_ID)
    #dir.create(str_c("../raw_ase/", subject_id), showWarnings = FALSE)
    #catch <- group_by(sample_read, TISSUE_ID) %>% group_map(SplitTissues, subject_id)
#}

#ase_files <- list.files("../wasp_counts", full.names = TRUE)

#map(ase_files, MapSplitTissues)

variants_keep <- c("3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant", "missense_variant", "synonymous_variant", "non_coding_transcript_exon_variant", "stop_gained")

#tissue_files <- list.files("../ase_by_tissue", full.names = TRUE, recursive = TRUE)
#tissue_names <- basename(tissue_files) %>% unique

tissue_files_reduce <- list.files("../ase_by_tissue_reduce", full.names = TRUE, recursive = TRUE) %>%
    str_subset("BLDDER|CVSEND|CVXECT|FLLPNT|KDNMDL", negate = TRUE)
tissue_reduce <- map(tissue_files_reduce, read_tsv)
names(tissue_reduce) <- basename(tissue_files_reduce)

#tissue_test <- read_tsv(tissue_files_reduce[[1]])
#tissue_test2 <- read_tsv("../ase_by_tissue/ADPSBQ")

AnnotateVariants <- function(tissue_table, tissue_name, variants_annot) {
    #tissue_table <- read_tsv(tissue_file)

    tissue_join <- inner_join(tissue_table, variants_annot)
    gc()
    write_tsv(tissue_join, str_c("../ase_by_tissue_annot/", tissue_name))
}

annot_filter <- filter(all_annot_reduce, is_in(Consequence, variants_keep)) 
temp <- AnnotateVariants(tissue_reduce[1], annot_filter)
imap(tissue_reduce, AnnotateVariants, annot_filter)

#old_annot <- read_tsv("../ase_by_tissue_annot_recent/ADPSBQ")
#new_annot <- read_tsv("../ase_by_tissue_annot/ADPSBQ")

tissue_files_annot <- list.files("../ase_by_tissue_annot", full.names = TRUE, recursive = TRUE) %>%
    str_subset("BLDDER|CVSEND|CVXECT|FLLPNT|KDNMDL", negate = TRUE)
tissue_annot <- map(tissue_files_annot, read_tsv)
names(tissue_annot) <- basename(tissue_files_annot)
saveRDS(tissue_annot, "tissue_annot.rda", compress = FALSE)

tissue_annot <- readRDS("./tissue_annot.rda")
#test_read <- read_tsv("../ase_by_tissue_annot_old/ADPSBQ")

GetProteinCoding <- function(tissue_annot, tissue_name) {
    tissue_stopgain <- filter(tissue_annot, Consequence == "stop_gained")
    tissue_stopgain
}

MergeASEPrediction <- function(tissue_annot, tissue_name, ase_predictions) {
    tissue_stopgain_prediction <- left_join(tissue_annot, ase_predictions) %>%
        filter(!is.na(Predicted_ASE_cat))
    tissue_stopgain_prediction
}

tissue_stopgain <- imap(tissue_annot, GetProteinCoding) %>% select

ase_predictions <- read_tsv("../../predictive_models/wasp_ase_predictions.txt") %>% 
    select(CHR, POS, REF_ALLELE, ALT_ALLELE, Predicted_ASE_cat) %>% distinct

tissue_stopgain_prediction <- imap(tissue_stopgain, MergeASEPrediction, ase_predictions)  %>% map(select, -Variant_ID, -TOTAL_COUNT, -Consequence)
map(tissue_stopgain_prediction, group_by, Predicted_ASE_cat) %>% map(tally)

GetCor <- function(variant_rows, variant_type) {
    ase_cor <- corAndPvalue(variant_rows$REF_RATIO1, variant_rows$REF_RATIO2)  %>%
        unlist %>% t %>% as_tibble
    ase_cor$Variant <- unlist(variant_type)
    ase_cor
}

tissue_pairs <- names(tissue_annot) %>% combn(2, simplify = FALSE)

TissuePair <- function(tissue_pair, tissue_list, variable_name) {
    tissue_pair1 <- select(tissue_list[[tissue_pair[1]]], -Variant_ID) 
    tissue_pair2 <- select(tissue_list[[tissue_pair[2]]], -Variant_ID) 
    
    colnames(tissue_pair1)[6] <- "TOTAL_COUNT1"
    colnames(tissue_pair2)[6] <- "TOTAL_COUNT2"
    colnames(tissue_pair1)[7] <- "REF_RATIO1"
    colnames(tissue_pair2)[7] <- "REF_RATIO2"

    two_tissue_join <- inner_join(tissue_pair1, tissue_pair2)
    two_tissue_filename <- str_c("tissue_join/", tissue_pair[1], "_", tissue_pair[2], ".rda")
    write_rds(two_tissue_join, two_tissue_filename)

    #if (nrow(two_tissue_join) > 0) {
        #tissue_pair_cor <- group_by(two_tissue_join, across(.cols = variable_name)) %>% group_map(GetCor) %>% bind_rows
        #tissue_pair_cor$Comb <- str_c(tissue_pair[1], "_vs_", tissue_pair[2])
    #} else {
        #num_variant_types <- unique(tissue_pair1[[variable_name]]) %>% length
        #print(num_variant_types)
        #tissue_pair_cor <- tibble(cor = rep(NA, num_variant_types),
            #p = rep(NA, num_variant_types),
            #Z = rep(NA, num_variant_types),
            #t = rep(NA, num_variant_types),
            #nObs = rep(NA, num_variant_types),
            #Variant = unique(tissue_pair1[[variable_name]]),
            #Comb = str_c(tissue_pair[1], "_vs_", tissue_pair[2]))
    #}
    #return(tissue_pair_cor)
}

tissue_cors <- map(tissue_pairs, TissuePair, tissue_annot, "Consequence") 
saveRDS(tissue_cors, "./tissue_cors_list.rda")
tissue_cors <- read_rds("./tissue_cors_list.rda")

tissue_cors_df <- bind_rows(tissue_cors) 
saveRDS(tissue_cors_df, "./tissue_cors_df_new.rda")

tissue_cors_nmd <- map(tissue_pairs, TissuePair, tissue_stopgain_prediction, "Predicted_ASE_cat")
tissue_cors_nmd_df <- bind_rows(tissue_cors_nmd)
tissue_cors_nmd_df$Variant %<>% str_replace_all("imbalanced", "predicted_NMD") %>% str_replace_all("balanced", "predicted_NMD_escape")

CollapseVariants <- function(df_list, bind_df) {
        df_read <- read_rds(df_list[1]) %>% select(CHR:ALT_ALLELE, Consequence)
        bind_df_distinct <- bind_rows(bind_df, df_read) %>% distinct

        if (length(df_list) > 1) {
            df_list_trim <- df_list[-1]
            CollapseVariants(df_list_trim, bind_df_distinct)
        } else {
            return(bind_df_distinct)
        }
}


tissue_pairs_list <- list.files("./tissue_join", full.names = T) 
read_df1 <- read_rds(tissue_pairs_list[1]) %>% select(CHR:ALT_ALLELE, Consequence)

combined_df <- CollapseVariants(tissue_pairs_list[2:length(tissue_pairs_list)], read_df1)
write_rds(combined_df, "./combined_df.rda")
combined_df <- read_rds("combined_df.rda")
combined_df$CHR %<>% str_remove_all("chr") %>% as.numeric

combined_df_pos <- select(combined_df, CHR, POS)
combined_df_pos$CHR %<>% str_remove_all("chr")
combined_df_pos$query <- str_c(combined_df_pos$CHR, ":", combined_df_pos$POS, "-", combined_df_pos$POS)
combined_df_out <- select(combined_df_pos, query)
write_tsv(combined_df_out, "ase_variants_query", col_names = FALSE)

gencode_annot <- readGFF("gencode.v26.annotation.gff3") %>% as.data.frame
gencode_cds_grange <- filter(gencode_annot, type == "CDS") %>%
    mutate(across(seqid, str_remove_all, "chr")) %>%
    select(seqid, start, end) %>% makeGRangesFromDataFrame

combined_df_cadd <- read_tsv("../../annotate_variants/ase_variants_cadd")
colnames(combined_df_cadd)[1:2] <- c("CHR", "POS")

combined_df_annot <- left_join(combined_df, combined_df_cadd) %>% distinct %>% 
    pivot_longer(GC:GerpS, names_to = "Variable", values_to = "Score") %>%
    filter(!(Variable == "mamPhyloP" & Score < -5)) %>%
    filter(!(Variable == "priPhyloP" & Score < -2.5)) %>%
    filter(!(Variable == "verPhyloP" & Score < -10))

combined_df_nc <- filter(combined_df_annot, Consequence == "non_coding_transcript_exon_variant") %>%
    select(CHR, POS) %>% distinct
colnames(combined_df_nc) <- c("seqid", "start")
combined_df_nc$end <- combined_df_nc$start
combined_df_nc_granges <- makeGRangesFromDataFrame(combined_df_nc)

nc_cds_overlaps <- findOverlaps(combined_df_nc_granges, gencode_cds_grange)

cadd_plot <- ggplot(combined_df_annot, aes(Consequence, Score)) +
    facet_wrap(~ Variable, ncol = 3, scales = "free_y") +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)
    )

CairoPDF("cadd_plot", width = 12, height = 10)
    print(cadd_plot)
dev.off()

#tissue_sex <- combn(c("PRSTTE", "TESTIS", "VAGINA", "UTERUS", "OVARY"), 2, simplify = FALSE) %>% reduce(rbind)
#tissue_sex_all <- rbind(tissue_sex, tissue_sex[,c(2,1)])

#tissue_cors_sex <- tibble(cor = NA, p = NA, Z = NA, t = NA, nObs = NA, Variant = NA, Tissue1 = tissue_sex_all[,1], Tissue2 = tissue_sex_all[,2]) %>%
    #filter(!(Tissue1 == "OVARY" & Tissue2 == "VAGINA" |
        #Tissue1 == "UTERUS" & Tissue2 == "VAGINA" |
        #Tissue1 == "OVARY" & Tissue2 == "UTERUS" |
        #Tissue1 == "PRSTTE" & Tissue2 == "TESTIS" |
        #Tissue2 == "OVARY" & Tissue1 == "VAGINA" |
        #Tissue2 == "UTERUS" & Tissue1 == "VAGINA" |
        #Tissue2 == "OVARU" & Tissue1 == "UTERUS" |
        #Tissue2 == "PRSTTE" & Tissue1 == "TESTIS"))

tissue_cors_df <- readRDS("./tissue_cors_df_new.rda") 
tissue_cors_all_df <- bind_rows(tissue_cors_df, tissue_cors_nmd_df) %>% arrange(Comb, Variant)

tissue_cors_all_df %<>% separate(Comb, into = c("Tissue1", "sep", "Tissue2"), sep = "_") %>% select(-sep)
tissue_cors_all_df$cor %<>% signif(digits = 3)
saveRDS(tissue_cors_all_df, "./tissue_cors_all_df.rda")

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

tissue_cors_5_prime_UTR <- filter(tissue_cors_df, Variant == "5_prime_UTR_variant")
colnames(tissue_cors_5_prime_UTR)[c(7,8)] <- c("Tissue2", "Tissue1")
tissue_cors_5_prime_UTR %<>% select(cor:Variant, Tissue1, Tissue2)

tissue_cors_stop_vs_5_prime_UTR <- bind_rows(tissue_cors_stop_gained, tissue_cors_5_prime_UTR, tissue_cors_self)
tissue_cors_sm_plot <- ggplot(tissue_cors_stop_vs_5_prime_UTR, aes(Tissue1, Tissue2, fill = cor, label = cor)) +
    geom_tile() +
    geom_text() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_fill_continuous(low = "white", high = muted("green"), na.value = "black", limits = c(0,1))

CairoPDF("cor_stop_vs_5_prime_UTR", height = 30, width = 30)
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

#non coding sanity checks
GetASEQuantiles <- function(ase_rows, consequence) {
    ase_quantiles <- quantile(ase_rows$REF_RATIO)
    ase_quantiles_df <- as_tibble(t(ase_quantiles)) %>%
        mutate(REF_RATIO_MAD = mad(ase_rows$REF_RATIO)) %>%
        mutate(MAD = mad(ase_rows$REF_RATIO)) %>%
        mutate(Consequence = unlist(consequence)) 
    ase_quantiles_df
}

MapGetASEQuantiles <- function(tissue_file) {
    tissue_file_df <- group_by(tissue_file, Consequence) %>% group_map(GetASEQuantiles) %>% bind_rows
    tissue_file_df
}

tissue_annot <- map(tissue_files_annot, read_tsv)
names(tissue_annot) <- basename(tissue_files_annot)


ase_quantiles_df <- bind_rows(ase_quantiles_list, .id = "Tissue")

MapGetCoverage <- function(tissue_table) {
    tissue_file_df <- group_by(tissue_table, Consequence) %>% summarise(Median = median(TOTAL_COUNT))
    tissue_file_df
}

ase_counts <- map(tissue_annot, MapGetCoverage) %>% bind_rows(.id = "Tissue")

tissue_annot_all <- map(tissue_annot, select, TOTAL_COUNT, REF_RATIO, Consequence) %>% bind_rows(.id = "Tissue")
write_rds(tissue_annot_all, "tissue_annot_all.rds")
tissue_annot_all <- read_rds("tissue_annot_all.rds")

ref_ratio_plot <- ggplot(tissue_annot_all, aes(Consequence, REF_RATIO, fill = Consequence)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ Tissue, ncol = 7, scales = "free_x") +
    theme_classic() +
    ylab("Reference ratio") +
    theme(panel.border = element_rect(color = "black", fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 30),
          axis.text.y = element_text(size = 27),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 27),
          legend.title = element_text(size = 30),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
    ) 


CairoPDF("refratio_boxplots", width = 42, height = 35)
    print(ref_ratio_plot)
dev.off()

tissue_annot_filter <- filter(tissue_annot_all, TOTAL_COUNT < 150)
counts_plot <- ggplot(tissue_annot_filter, aes(Consequence, TOTAL_COUNT, fill = Consequence)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~ Tissue, ncol = 7, scales = "free_x") +
    theme_classic() +
    ylab("Total counts") +
    theme(panel.border = element_rect(color = "black", fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 30),
          axis.text.y = element_text(size = 27),
          axis.title.y = element_text(size = 30),
          legend.text = element_text(size = 27),
          legend.title = element_text(size = 30),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
    )

CairoPDF("counts_boxplots", width = 42, height = 35)
    print(counts_plot)
dev.off()

objects_size <- map(ls(), function(thing) print(object.size(get(thing)), units = 'auto')) 
names(objects_size) <- ls()
unlist(objects_size) %>% sort
