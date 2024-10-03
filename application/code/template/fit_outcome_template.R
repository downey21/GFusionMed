
# -*- coding: utf-8 -*-

rm(list = ls())

pp = "pathway_information"

drug_name = "drug_name"

suppressPackageStartupMessages({
    library(GFusionMed)
})

source("/root/Project/GFusionMed/application/code/functions_for_data_analysis.R")

path_data <- "/root/Project/GFusionMed/application/data"
path_result <- "/root/Project/GFusionMed/application/result"

data_for_outcome <-
    load_data_for_outcome(
        path_data = path_data,
        pp = pp,
        drug_name = drug_name
    )

result_fit <- GFusionMed::fit_outcome_model(data_for_outcome)

save(result_fit, file = paste(path_result, "/fit_outcome_", drug_name, "_pp", pp, ".RData", sep = ""))

cat("End")
