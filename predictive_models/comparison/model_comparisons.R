library(modEvA)
library(caret)
library(ROCR)
library(pROC)

library(Cairo)
library(magrittr)
library(tidyverse)

#Lindeboom
lindeboom_modase <- readRDS("../lindeboom/glmnet_fit_modase.rda")
lindeboom_modase_df <- tibble(Lambda = lindeboom_modase$lambda, CVM = lindeboom_modase$cvm, CVSD = lindeboom_modase$cvsd, Model = "Lindeboom", Outcome = "MODASE")

lindeboom_refratio <- readRDS("../lindeboom/glmnet_fit_refratio.rda")
lindeboom_refratio_df <- tibble(Lambda = lindeboom_refratio$lambda, CVM = lindeboom_refratio$cvm, CVSD = lindeboom_refratio$cvsd, Model = "Lindeboom", Outcome = "REF_RATIO")

lindeboom_modase_cont <- readRDS("../lindeboom/glmnet_fit_modase_cont.rda")
lindeboom_modase_cont_df <- tibble(Lambda = lindeboom_modase_cont$lambda, CVM = lindeboom_modase_cont$cvm, CVSD = lindeboom_modase_cont$cvsd, Model = "Lindeboom", Outcome = "MODASE")

lindeboom_refratio_cont <- readRDS("../lindeboom/glmnet_fit_refratio_cont.rda")
lindeboom_refratio_cont_df <- tibble(Lambda = lindeboom_refratio_cont$lambda, CVM = lindeboom_refratio_cont$cvm, CVSD = lindeboom_refratio_cont$cvsd, Model = "Lindeboom", Outcome = "REF_RATIO")

#Full model
full_modase <- readRDS("../full/glmnet_fit_modase.rda")
full_modase_df <- tibble(Lambda = full_modase$lambda, CVM = full_modase$cvm, CVSD = full_modase$cvsd, Model = "No Tissue", Outcome = "MODASE")

full_refratio <- readRDS("../full/glmnet_fit_refratio.rda")
full_refratio_df <- tibble(Lambda = full_refratio$lambda, CVM = full_refratio$cvm, CVSD = full_refratio$cvsd, Model = "No Tissue", Outcome = "REF_RATIO")

full_modase_cont <- readRDS("../full/glmnet_fit_modase_cont.rda")
full_modase_cont_df <- tibble(Lambda = full_modase_cont$lambda, CVM = full_modase_cont$cvm, CVSD = full_modase_cont$cvsd, Model = "No Tissue", Outcome = "MODASE")

full_refratio_cont <- readRDS("../full/glmnet_fit_refratio_cont.rda")
full_refratio_cont_df <- tibble(Lambda = full_refratio_cont$lambda, CVM = full_refratio_cont$cvm, CVSD = full_refratio_cont$cvsd, Model = "No Tissue", Outcome = "REF_RATIO")

#Tissue
tissue_refratio <- readRDS("../tissue/glmnet_fit_refratio.rda")
tissue_refratio_df <- tibble(Lambda = tissue_refratio$lambda, CVM = tissue_refratio$cvm, CVSD = tissue_refratio$cvsd, Model = "Tissue", Outcome = "REF_RATIO")

tissue_refratio_cont <- readRDS("../tissue/glmnet_fit_refratio_cont.rda")
tissue_refratio_cont_df <- tibble(Lambda = tissue_refratio_cont$lambda, CVM = tissue_refratio_cont$cvm, CVSD = tissue_refratio_cont$cvsd, Model = "Tissue", Outcome = "REF_RATIO")

#Full model
subject_modase <- readRDS("../subject/glmnet_fit_modase.rda")
subject_modase_df <- tibble(Lambda = subject_modase$lambda, CVM = subject_modase$cvm, CVSD = subject_modase$cvsd, Model = "Subject", Outcome = "MODASE")

subject_refratio <- readRDS("../subject/glmnet_fit_refratio.rda")
subject_refratio_df <- tibble(Lambda = subject_refratio$lambda, CVM = subject_refratio$cvm, CVSD = subject_refratio$cvsd, Model = "Subject", Outcome = "REF_RATIO")

subject_modase_cont <- readRDS("../subject/glmnet_fit_modase_cont.rda")
subject_modase_cont_df <- tibble(Lambda = subject_modase_cont$lambda, CVM = subject_modase_cont$cvm, CVSD = subject_modase_cont$cvsd, Model = "Subject", Outcome = "MODASE")

subject_refratio_cont <- readRDS("../subject/glmnet_fit_refratio_cont.rda")
subject_refratio_cont_df <- tibble(Lambda = subject_refratio_cont$lambda, CVM = subject_refratio_cont$cvm, CVSD = subject_refratio_cont$cvsd, Model = "Subject", Outcome = "REF_RATIO")

#Combined
refratio_combined <- bind_rows(lindeboom_refratio_df, full_refratio_df, tissue_refratio_df, subject_refratio_df)
saveRDS(refratio_combined, "./refratio_combined.rda")
refratio_cont_combined <- bind_rows(lindeboom_refratio_cont_df, full_refratio_cont_df, tissue_refratio_cont_df, subject_refratio_cont_df)
saveRDS(refratio_cont_combined, "./refratio_cont_combined.rda")

modase_combined <- bind_rows(lindeboom_modase_df, full_modase_df, subject_modase_df)
saveRDS(modase_combined, "./modase_combined.rda")
modase_cont_combined <- bind_rows(lindeboom_modase_cont_df, full_modase_cont_df, subject_modase_cont_df)
saveRDS(modase_cont_combined, "./modase_cont_combined.rda")

cat_combined <- bind_rows(refratio_combined, modase_combined) %>% filter(Lambda < 0.15)
cont_combined <- bind_rows(refratio_cont_combined, modase_cont_combined) 

PlotPerformance_Simple <- function(glmnet_cvm_combined_df, filename, ylabel) {
    glmnet_cvm_combined_df$Ymin <- glmnet_cvm_combined_df$CVM - glmnet_cvm_combined_df$CVSD
    glmnet_cvm_combined_df$Ymax <- glmnet_cvm_combined_df$CVM + glmnet_cvm_combined_df$CVSD

    cvm_plot <- ggplot(glmnet_cvm_combined_df, aes(log10(Lambda), CVM, color = Model, ymin = Ymin, ymax = Ymax)) +
        geom_line() +
        geom_errorbar(width = 0) +
        theme_classic() +
        scale_x_reverse() +
        ylab(ylabel) +
        xlab(expression(paste(log[10], " Lambda"))) +
        theme(plot.background = element_blank(),
              panel.border = element_rect(fill = NA),
            legend.background = element_blank())

    CairoPDF(filename, height = 4, width = 7)
        print(cvm_plot)
    dev.off()
}

PlotPerformance <- function(glmnet_cvm_combined_df, filename, ylabel) {
    glmnet_cvm_combined_df$Ymin <- glmnet_cvm_combined_df$CVM - glmnet_cvm_combined_df$CVSD
    glmnet_cvm_combined_df$Ymax <- glmnet_cvm_combined_df$CVM + glmnet_cvm_combined_df$CVSD

    cvm_plot <- ggplot(glmnet_cvm_combined_df, aes(log10(Lambda), CVM, color = Model, ymin = Ymin, ymax = Ymax)) +
        geom_line() +
        geom_errorbar(width = 0) +
        facet_wrap(~ Outcome, ncol = 2) +
        theme_classic() +
        scale_x_reverse() +
        ylab(ylabel) +
        xlab(expression(paste(log[10], " Lambda"))) +
        theme(plot.background = element_blank(),
              panel.border = element_rect(fill = NA),
            legend.background = element_blank())

    CairoPDF(filename, height = 4, width = 14)
        print(cvm_plot)
    dev.off()
}

#PlotPerformance(refratio_combined, "refratio_cvm")
#PlotPerformance(refratio_cont_combined, "refratio_cont_cvm")

#PlotPerformance(modase_combined, "modase_cvm")
#PlotPerformance(modase_cont_combined, "modase_cont_cvm")

PlotPerformance(cat_combined, "cat_cvm", "AUC")
PlotPerformance(cont_combined, "cont_cvm", "MSE")

refratio_combined_filter <- filter(refratio_combined, Model != "Subject") %>%
    filter(Lambda > 0.003)
PlotPerformance_Simple(refratio_combined_filter, "refratio_cvm_simple", "AUC")

