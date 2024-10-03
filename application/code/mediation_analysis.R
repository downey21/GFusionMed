
# -*- coding: utf-8 -*-

rm(list = ls())

suppressPackageStartupMessages({
    library(GFusionMed)
})

source("/root/Project/GFusionMed/application/code/functions_for_data_analysis.R")

path_data <- "/root/Project/GFusionMed/application/data"
path_result <- "/root/Project/GFusionMed/application/result"
path_result_analysis <- "/root/Project/GFusionMed/application/result_analysis"

vector_drug <- c("erlotinib", "gefitinib", "afatinib", "dacomitinib", "osimertinib", "cisplatin")

result <- mediation_summary_all(path_data, path_result, vector_drug)

write.csv(result, file = paste0(path_result_analysis, "/mediation_analysis_results.csv"), row.names = FALSE)
