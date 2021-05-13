library(rtracklayer)
library(magrittr)
library(tidyverse)

gtf_in <- readGFF("./gencode.v26.annotation.gtf")
gtf_in$Gene <- str_remove_all(gtf_in$gene_id, "\\..*$")

gtf_exons <- filter(gtf_in, type == "exon")
gtf_exons_select <- select(gtf_exons, seqid, start, end, transcript_id, exon_number) %>% as_tibble
colnames(gtf_exons_select)[1:5] <- c("Chrom", "exon_start", "exon_end", "Feature", "Exon_position")
gtf_exons_select$Feature %<>% str_remove_all("\\..*$")
gtf_exons_select$Chrom %<>% as.integer
gtf_exons_select$Exon_position %<>% as.integer

variants_load <- read_tsv("./vep_stopgains_annot.txt", skip = 89, col_types = str_c(rep("c", 67), collapse = "")) %>% filter(grepl("stop_gained", Consequence) & CANONICAL == "YES")

variants_exon_split <- str_split_fixed(variants_load$EXON, "\\/", 2)
variants_load$Exon_position <- as.integer(variants_exon_split[,1])
variants_load$N_Exons <- as.integer(variants_exon_split[,2])

variants_cdna_split <- str_split_fixed(variants_load$cDNA_position, "/", 2)
variants_cds_split <- str_split_fixed(variants_load$CDS_position, "/", 2)
variants_protein_split <- str_split_fixed(variants_load$Protein_position, "/", 2)

variants_load$cDNA_position_only <- as.integer(variants_cdna_split[,1])
variants_load$cDNA_length <- as.integer(variants_cdna_split[,2])
variants_load$relative_cDNA_position <- variants_load$cDNA_position_only / variants_load$cDNA_length

variants_load$CDS_position_only <- as.integer(variants_cds_split[,1])
variants_load$CDS_length <- as.integer(variants_cds_split[,2])
variants_load$relative_CDS_position <- variants_load$CDS_position_only / variants_load$CDS_length

variants_load$Protein_position_only <- as.integer(variants_protein_split[,1])
variants_load$Protein_length <- as.integer(variants_protein_split[,2])
variants_load$relative_Protein_position <- variants_load$Protein_position_only / variants_load$Protein_length

variants_lofinfo_split <- str_split_fixed(variants_load$LoF_info, ",", 5) %>% as_tibble %>%
    mutate_all(str_remove_all, "^.*:") %>%
    mutate_at(vars(V1:V4), as.numeric)
colnames(variants_lofinfo_split) <- c("PERCENTILE", "GERP_DIST", "BP_DIST", "DIST_FROM_LAST_EXON", "50_BP_RULE")

variants_load_lof <- cbind(variants_load, variants_lofinfo_split) %>% as_tibble

variants_load_97 <- read_tsv("./vep_stopgains_annot_97.txt", skip = 51, col_types = str_c(rep("c", 29), collapse = ""))
variants_load_gnomad <- select(variants_load_97, Location, Allele, Gene, gnomAD_AF) %>% distinct 

variants_load_join <- left_join(variants_load_lof, variant_load_gnomad) 
variants_load_join$gnomAD_AF[variants_load_join$gnomAD_AF == "-"] <- "0"
variants_load_join$gnomAD_AF[is.na(variants_load_join$gnomAD_AF)] <- 0
variants_load_join$gnomAD_AF %<>% as.numeric

variants_location_split <- str_split_fixed(variants_load_join$Location, ":", 2)
variants_load_join$Chrom <- as.integer(variants_location_split[,1])
variants_load_join$Pos <- as.integer(variants_location_split[,2])
colnames(variants_load_join)[3] <- "Alt"

variants_load_join2 <- left_join(variants_load_join, gtf_exons_select)
variants_load_join2$Exon_Length <- variants_load_join2$exon_end - variants_load_join2$exon_start

cadd_variants <- read_tsv("./stopgain_variants_grepped_from_cadd.txt")
colnames(cadd_variants)[1] <- "Chrom"
colnames(cadd_variants)[19] <- "Gene"
cadd_select <- filter(cadd_variants, AnnoType == "CodingTranscript") %>% select(Chrom:Alt, Gene, relcDNApos, relCDSpos, relProtPos, minDistTSS, minDistTSE, GerpN, GerpS, GC, CpG, priPhCons:verPhyloP, PHRED, contains("sum"))
cadd_select$GerpN[is.na(cadd_select$GerpN)] <- 0
cadd_select$GerpS[is.na(cadd_select$GerpS)] <- 0

variants_load_cadd <- left_join(variants_load_join2, cadd_select)
saveRDS(variants_load_cadd, "./variants_annotated.rda")
