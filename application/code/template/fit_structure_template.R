
# -*- coding: utf-8 -*-

rm(list = ls())

pp = "pathway_information"

suppressPackageStartupMessages({
    library(GFusionMed)
})

source("/root/Project/GFusionMed/application/code/functions_for_data_analysis.R")

path_data <- "/root/Project/GFusionMed/application/data"
path_result <- "/root/Project/GFusionMed/application/result"

data_for_structure <-
    load_data_for_structure(
        path_data = path_data,
        pp = pp
    )

result_fit <- GFusionMed::fit_structure_model(data_for_structure, cores = 3)

save(result_fit, file = paste(path_result, "/fit_structure_pp", pp, ".RData", sep = ""))

cat("End")

# example_data_for_structure <- data_for_structure
# save(example_data_for_structure, file = file.path(path_result, "example_data_for_structure.rda"), version = 2)
