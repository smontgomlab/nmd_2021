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
#DTHHRDY dropped due to missingness
glm_fit_data_refratio <- select(wasp_final_noagg, CHR, REF_RATIO_cat, TOTAL_COUNT, `50_BP_RULE`, Dist_to_End, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
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
glmnet_fit_min <- which.max(glmnet_fit_refratio$cvm)
glmnet_fit_coefs <- coef(glmnet_fit_refratio$glmnet.fit)[,glmnet_fit_min]

saveRDS(glmnet_fit_refratio, "glmnet_fit_refratio.rda")

#MODASE_cat
glm_fit_data_modase <- select(wasp_final_noagg, CHR, MODASE_cat, TOTAL_COUNT, `50_BP_RULE`, Dist_to_End, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
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
glmnet_fit_min_modase <- which.max(glmnet_fit_modase$cvm)
glmnet_fit_coefs_modase <- coef(glmnet_fit_modase$glmnet.fit)[,glmnet_fit_min_modase]
saveRDS(glmnet_fit_modase, "glmnet_fit_modase.rda")
glmnet_fit_modase_predict_prob <- predict(glmnet_fit_modase$glmnet.fit, glmnet_matrix_modase, type = "response")[,glmnet_fit_min_modase]
glmnet_fit_modase_predict_class <- predict(glmnet_fit_modase$glmnet.fit, glmnet_matrix_modase, type = "class")[,glmnet_fit_min_modase]
ase_predictions <- filter(wasp_final_noagg, !is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY)) %>% 
    mutate(Predicted_MODASE_cat_prob = glmnet_fit_modase_predict_prob, Predicted_MODASE_cat = glmnet_fit_modase_predict_class)
write_tsv(ase_predictions, "ase_predictions.txt")

glmnet_fixedinf <- fixedLassoInf(glmnet_matrix_modase, as.integer(factor(glm_fit_data_modase$MODASE_cat)), coef(glmnet_fit_modase$glmnet.fit, glmnet_fit_modase$lambda.min), glmnet_fit_modase$lambda.min, family = "binomial")
glmnet_fixedinf_df <- tibble(Variable = names(glmnet_fixedinf$vars), 
    Coefficient = formatC(exp(glmnet_fixedinf$coef0), format = "f", digits = 2),
    Coefficient_numeric = exp(glmnet_fixedinf$coef0),
    CI_LB = formatC(exp(glmnet_fixedinf$ci[,1]), format = "f", digits = 2),
    CI_UB = formatC(exp(glmnet_fixedinf$ci[,2]), format = "f", digits = 2),
    Pvalue = glmnet_fixedinf$pv,
    Zscore = glmnet_fixedinf$zscore0) %>%
    mutate_if(is.numeric, signif, digits = 3) %>% arrange(Pvalue)

saveRDS(glmnet_fixedinf_df, "./glmnet_inference_df.rda")

#REF_RATIO_cat
glm_fit_data_refratio_cont <- select(wasp_final_noagg, CHR, logit_REF_RATIO, TOTAL_COUNT, `50_BP_RULE`, Dist_to_End, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

glm_fit_data_refratio_cont$`50_BP_RULE` <- glm_fit_data_refratio_cont$`50_BP_RULE` == "FAIL"

glm_fit_nochr_refratio_cont <- select(glm_fit_data_refratio_cont, -CHR, -logit_REF_RATIO)
glmnet_matrix_refratio_cont <- scale(model.matrix(~ ., glm_fit_nochr_refratio_cont)[,-1], center = TRUE)
glmnet_fit_refratio_cont <- cv.glmnet(glmnet_matrix_refratio_cont, scale(glm_fit_data_refratio_cont$logit_REF_RATIO), foldid = as.integer(factor(glm_fit_data_refratio_cont$CHR)), family = "gaussian", alpha = 1.0, type.measure = "mse")
glmnet_fit_min_refratio_cont <- which.min(glmnet_fit_refratio_cont$cvm)
glmnet_fit_coefs_refratio_cont <- coef(glmnet_fit_refratio_cont$glmnet.fit)[,glmnet_fit_min_refratio_cont]

saveRDS(glmnet_fit_refratio_cont, "glmnet_fit_refratio_cont.rda")

#MODASE_cat
glm_fit_data_modase_cont <- select(wasp_final_noagg, CHR, logit_MODASE, TOTAL_COUNT, `50_BP_RULE`, Dist_to_End, Long_Exon, Near_Start, LoF, cDNA_position_only, CDS_position_only, Protein_position_only, GerpN:PHRED, gnomAD_AF, SMRIN, SMTSISCH, SEX, AGE, DTHHRDY) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) %>%
    filter(!is.na(GC)) %>%
    filter(!is.na(SMTSISCH)) %>%
    filter(!is.na(DTHHRDY))

glm_fit_data_modase_cont$`50_BP_RULE` <- glm_fit_data_modase_cont$`50_BP_RULE` == "FAIL"

glm_fit_nochr_modase_cont <- select(glm_fit_data_modase_cont, -CHR, -logit_MODASE)
glmnet_matrix_modase_cont <- scale(model.matrix(~ ., glm_fit_nochr_modase_cont)[,-1], center = TRUE)
glmnet_fit_modase_cont <- cv.glmnet(glmnet_matrix_modase_cont, scale(glm_fit_data_modase_cont$logit_MODASE), foldid = as.integer(factor(glm_fit_data_modase_cont$CHR)), family = "gaussian", alpha = 1.0, type.measure = "mse")
glmnet_fit_min_modase_cont <- which.min(glmnet_fit_modase_cont$cvm)
glmnet_fit_coefs_modase_cont <- coef(glmnet_fit_modase_cont$glmnet.fit)[,glmnet_fit_min_modase_cont]
saveRDS(glmnet_fit_modase_cont, "glmnet_fit_modase_cont.rda")
