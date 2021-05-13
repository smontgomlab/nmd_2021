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

#REF_RATIO_cat
glm_fit_data_refratio <- select(wasp_final_noagg, CHR, REF_RATIO_cat, `50_BP_RULE`, Long_Exon, Near_Start) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) 

glm_fit_data_refratio$`50_BP_RULE` <- glm_fit_data_refratio$`50_BP_RULE` == "FAIL"
glm_fit_data_refratio$REF_RATIO_cat %<>% factor(levels = c("balanced", "imbalanced"))

glm_fit_nochr_refratio <- select(glm_fit_data_refratio, -CHR, -REF_RATIO_cat)
glmnet_matrix_refratio <- scale(model.matrix(~ ., glm_fit_nochr_refratio)[,-1], center = TRUE)
glmnet_fit_refratio <- cv.glmnet(glmnet_matrix_refratio, glm_fit_data_refratio$REF_RATIO_cat, foldid = as.integer(factor(glm_fit_data_refratio$CHR)), family = "binomial", alpha = 1.0, type.measure = "auc")
saveRDS(glmnet_fit_refratio, "glmnet_fit_refratio.rda")

best_error_refratio <- which.max(glmnet_fit_refratio$cvm)
glmnet_pred_refratio <- predict(glmnet_fit_refratio$glmnet.fit, glmnet_matrix_refratio, type = "response")
glmnet_pred_refratio_best <- glmnet_pred_refratio[,best_error_refratio]

wasp_nomissing <- filter(wasp_final_noagg, !is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) 
wasp_nomissing$Predicted_ASE_lindeboom <- glmnet_pred_refratio_best
wasp_nomissing$Predicted_ASE_lindeboom_cat <- round(wasp_nomissing$Predicted_ASE_lindeboom) %>%
    factor %>%
    fct_recode(imbalanced = "1", balanced = "0") %>%
    factor(levels = c("balanced", "imbalanced"))   

write_tsv(wasp_nomissing, "./wasp_ase_predictions_lindeboom_refratio.txt")

should_imbalanced <- filter(wasp_nomissing, `50_BP_RULE` == "PASS" & Long_Exon == FALSE & CDS_position_only >= 150)
improper_balanced <- filter(should_imbalanced, Predicted_ASE_lindeboom_cat == "balanced") 
proper_imbalanced <- filter(should_imbalanced, Predicted_ASE_lindeboom_cat == "imbalanced") 

test_variant <- filter(wasp_nomissing, SYMBOL == "HBB")

#MODASE_cat
glm_fit_data_modase <- select(wasp_final_noagg, CHR, MODASE_cat, `50_BP_RULE`, Long_Exon, Near_Start) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) 

glm_fit_data_modase$`50_BP_RULE` <- glm_fit_data_modase$`50_BP_RULE` == "FAIL"
glm_fit_data_modase$MODASE_cat %<>% factor(levels = c("balanced", "imbalanced"))

glm_fit_nochr_modase <- select(glm_fit_data_modase, -CHR, -MODASE_cat)
glmnet_matrix_modase <- scale(model.matrix(~ ., glm_fit_nochr_modase)[,-1], center = TRUE)
glmnet_fit_modase <- cv.glmnet(glmnet_matrix_modase, glm_fit_data_modase$MODASE_cat, foldid = as.integer(factor(glm_fit_data_modase$CHR)), family = "binomial", alpha = 1.0, type.measure = "auc")
saveRDS(glmnet_fit_modase, "glmnet_fit_modase.rda")

#REF_RATIO_cat
glm_fit_data_refratio_cont <- select(wasp_final_noagg, CHR, logit_REF_RATIO, `50_BP_RULE`, Long_Exon, Near_Start) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) 

glm_fit_data_refratio_cont$`50_BP_RULE` <- glm_fit_data_refratio_cont$`50_BP_RULE` == "FAIL"

glm_fit_nochr_refratio_cont <- select(glm_fit_data_refratio_cont, -CHR, -logit_REF_RATIO)
glmnet_matrix_refratio_cont <- scale(model.matrix(~ ., glm_fit_nochr_refratio_cont)[,-1], center = TRUE)
glmnet_fit_refratio_cont <- cv.glmnet(glmnet_matrix_refratio_cont, scale(glm_fit_data_refratio_cont$logit_REF_RATIO), foldid = as.integer(factor(glm_fit_data_refratio_cont$CHR)), family = "gaussian", alpha = 1.0, type.measure = "mse")
saveRDS(glmnet_fit_refratio_cont, "glmnet_fit_refratio_cont.rda")

#MODASE_cat
glm_fit_data_modase_cont <- select(wasp_final_noagg, CHR, logit_MODASE, `50_BP_RULE`, Long_Exon, Near_Start) %>% 
    filter(!is.na(`50_BP_RULE`)) %>%
    filter(nchar(`50_BP_RULE`) > 0) 

glm_fit_data_modase_cont$`50_BP_RULE` <- glm_fit_data_modase_cont$`50_BP_RULE` == "FAIL"

glm_fit_nochr_modase_cont <- select(glm_fit_data_modase_cont, -CHR, -logit_MODASE)
glmnet_matrix_modase_cont <- scale(model.matrix(~ ., glm_fit_nochr_modase_cont)[,-1], center = TRUE)
glmnet_fit_modase_cont <- cv.glmnet(glmnet_matrix_modase_cont, scale(glm_fit_data_modase_cont$logit_MODASE), foldid = as.integer(factor(glm_fit_data_modase_cont$CHR)), family = "gaussian", alpha = 1.0, type.measure = "mse")
saveRDS(glmnet_fit_modase_cont, "glmnet_fit_modase_cont.rda")
