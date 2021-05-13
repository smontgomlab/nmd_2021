library(modEvA)
library(MASS)
library(caret)
library(glmnet)
library(ROCR)
library(pROC)
library(parallel)
library(selectiveInference)
library(limma)
library(lme4)
library(contrast)

library(broom)
library(forestplot)
library(rlist)
library(magrittr)
library(tidyverse)

wasp_final_noagg <- readRDS("../wasp_final_noagg.rda") 
wasp_final_noagg$Dist_to_End <- wasp_final_noagg$CDS_length - wasp_final_noagg$CDS_position_only

#REF_RATIO_cat
glm_fit_data_refratio <- select(wasp_final_noagg, CHR, REF_RATIO_cat, TOTAL_COUNT, `50_BP_RULE`, Dist_to_End, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, TISSUE_ID, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

#glm_fit_data_refratio$`50_BP_RULE` <- glm_fit_data_refratio$`50_BP_RULE` == "FAIL"
glm_fit_data_refratio$REF_RATIO_cat %<>% factor(levels = c("balanced", "imbalanced"))

glm_fit_data_nmd <- filter(glm_fit_data_refratio, `50_BP_RULE` == "PASS") %>% select(CHR, REF_RATIO_cat, TISSUE_ID)
glm_fit_data_nmd <- filter(glm_fit_data_refratio, `50_BP_RULE` == "PASS" & Long_Exon == FALSE & Near_Start == FALSE) %>% select(CHR, REF_RATIO_cat, TISSUE_ID)

glm_fit_nochr_refratio <- select(glm_fit_data_nmd, -CHR, -REF_RATIO_cat)
glmnet_matrix_refratio <- scale(model.matrix(~ ., glm_fit_nochr_refratio)[,-1], center = TRUE)
glmnet_fit_refratio <- cv.glmnet(glmnet_matrix_refratio, glm_fit_data_nmd$REF_RATIO_cat, foldid = as.integer(factor(glm_fit_data_nmd$CHR)), family = "binomial", alpha = 1.0, type.measure = "auc")
saveRDS(glmnet_fit_refratio, "glmnet_fit_refratio.rda")

glmnet_cvm_combined_df <- tibble(Lambda = glmnet_fit_refratio$lambda, CVM = glmnet_fit_refratio$cvm, CVSD = glmnet_fit_refratio$cvsd)
glmnet_cvm_combined_df$Ymin <- glmnet_cvm_combined_df$CVM - glmnet_cvm_combined_df$CVSD
glmnet_cvm_combined_df$Ymax <- glmnet_cvm_combined_df$CVM + glmnet_cvm_combined_df$CVSD

cvm_plot <- ggplot(glmnet_cvm_combined_df, aes(log10(Lambda), CVM, ymin = Ymin, ymax = Ymax)) +
    geom_line() +
    geom_errorbar(width = 0) +
    theme_classic() +
    scale_x_reverse() +
    ylab("AUC") +
    xlab(expression(paste(log[10], " Lambda"))) +
    theme(plot.background = element_blank(),
          panel.border = element_rect(fill = NA),
        legend.background = element_blank())

CairoPDF("refratio_50nt_only", height = 4, width = 7)
    print(cvm_plot)
dev.off()
