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

subject_var_counts <- group_by(wasp_final_noagg, SUBJECT_ID, CHR) %>% 
    tally %>% group_by(SUBJECT_ID) %>% tally
subject_single_chrom <- filter(subject_var_counts, n == 1)

glm_fit_data_subject <- filter(wasp_final_noagg, !is_in(SUBJECT_ID, subject_single_chrom$SUBJECT_ID))

#REF_RATIO_cat
glm_fit_data_refratio <- select(glm_fit_data_subject, CHR, REF_RATIO_cat, TOTAL_COUNT, `50_BP_RULE`, Dist_to_End, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SUBJECT_ID, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

glm_fit_data_refratio$`50_BP_RULE` <- glm_fit_data_refratio$`50_BP_RULE` == "FAIL"
glm_fit_data_refratio$REF_RATIO_cat %<>% factor(levels = c("balanced", "imbalanced"))

glm_fit_nochr_refratio <- select(glm_fit_data_refratio, -CHR, -REF_RATIO_cat)
glmnet_matrix_refratio <- scale(model.matrix(~ ., glm_fit_nochr_refratio)[,-1], center = TRUE)
glmnet_fit_refratio <- cv.glmnet(glmnet_matrix_refratio, glm_fit_data_refratio$REF_RATIO_cat, foldid = as.integer(factor(glm_fit_data_refratio$CHR)), family = "binomial", alpha = 1.0, type.measure = "auc")
saveRDS(glmnet_fit_refratio, "glmnet_fit_refratio.rda")

#MODASE_cat
glm_fit_data_modase <- select(glm_fit_data_subject, CHR, MODASE_cat, TOTAL_COUNT, `50_BP_RULE`, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SUBJECT_ID, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

glm_fit_data_modase$`50_BP_RULE` <- glm_fit_data_modase$`50_BP_RULE` == "FAIL"
glm_fit_data_modase$MODASE_cat %<>% factor(levels = c("balanced", "imbalanced"))

glm_fit_nochr_modase <- select(glm_fit_data_modase, -CHR, -MODASE_cat)
glmnet_matrix_modase <- scale(model.matrix(~ ., glm_fit_nochr_modase)[,-1], center = TRUE)
glmnet_fit_modase <- cv.glmnet(glmnet_matrix_modase, glm_fit_data_modase$MODASE_cat, foldid = as.integer(factor(glm_fit_data_modase$CHR)), family = "binomial", alpha = 1.0, type.measure = "auc")
saveRDS(glmnet_fit_modase, "glmnet_fit_modase.rda")

#REF_RATIO_cat
glm_fit_data_refratio_cont <- select(glm_fit_data_subject, CHR, logit_REF_RATIO, TOTAL_COUNT, `50_BP_RULE`, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SUBJECT_ID, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

glm_fit_data_refratio_cont$`50_BP_RULE` <- glm_fit_data_refratio_cont$`50_BP_RULE` == "FAIL"

glm_fit_nochr_refratio_cont <- select(glm_fit_data_refratio_cont, -CHR, -logit_REF_RATIO)
glmnet_matrix_refratio_cont <- scale(model.matrix(~ ., glm_fit_nochr_refratio_cont)[,-1], center = TRUE)
glmnet_fit_refratio_cont <- cv.glmnet(glmnet_matrix_refratio_cont, scale(glm_fit_data_refratio_cont$logit_REF_RATIO), foldid = as.integer(factor(glm_fit_data_refratio_cont$CHR)), family = "gaussian", alpha = 1.0, type.measure = "mse")
saveRDS(glmnet_fit_refratio_cont, "glmnet_fit_refratio_cont.rda")

#MODASE_cat
glm_fit_data_modase_cont <- select(glm_fit_data_subject, CHR, logit_MODASE, TOTAL_COUNT, `50_BP_RULE`, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SUBJECT_ID, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

glm_fit_data_modase_cont$`50_BP_RULE` <- glm_fit_data_modase_cont$`50_BP_RULE` == "FAIL"

glm_fit_nochr_modase_cont <- select(glm_fit_data_modase_cont, -CHR, -logit_MODASE)
glmnet_matrix_modase_cont <- scale(model.matrix(~ ., glm_fit_nochr_modase_cont)[,-1], center = TRUE)
glmnet_fit_modase_cont <- cv.glmnet(glmnet_matrix_modase_cont, scale(glm_fit_data_modase_cont$logit_MODASE), foldid = as.integer(factor(glm_fit_data_modase_cont$CHR)), family = "gaussian", alpha = 1.0, type.measure = "mse")
saveRDS(glmnet_fit_modase_cont, "glmnet_fit_modase_cont.rda")
