
# -*- encoding: utf-8 -*-

rm(list = ls())

suppressPackageStartupMessages({
    library(GFusionMed)
})

path_result <- "/root/Project/GFusionMed/application/result"
path_result_analysis <- "/root/Project/GFusionMed/application/result_analysis"

pp <- 9 # RTK Pathway
drug_name <- "erlotinib"
exposure <- "mRNA_EGFR"

load(paste(path_result, "/fit_structure_pp", pp, ".RData", sep = ""))
result_structure <- result_fit

load(paste(path_result,"/fit_outcome_", drug_name, "_pp", pp, ".RData", sep = ""))
result_outcome <- result_fit

RTK_erlotinib_mRNA_EGFR <- GFusionMed::perform_mediation_analysis(result_structure, result_outcome, exposure)

str(RTK_erlotinib_mRNA_EGFR)

print(RTK_erlotinib_mRNA_EGFR)

GFusionMed::plot_network(
    result_structure = result_structure,
    path = path_result_analysis,
    file_name = "plot_network_structure_1"
)

GFusionMed::plot_network(
    result_outcome = result_outcome,
    path = path_result_analysis,
    file_name = "plot_network_outcome_1"
)

GFusionMed::plot_network(
    result_structure = result_structure,
    result_outcome = result_outcome,
    path = path_result_analysis,
    file_name = "plot_network_structure_outcome_1"
)

GFusionMed::plot_network(
    result_structure = result_structure,
    result_outcome = result_outcome,
    exposure = exposure,
    path = path_result_analysis,
    file_name = "plot_network_structure_outcome_exposure_1"
)

GFusionMed::plot_network(
    result_structure = result_structure,
    result_outcome = result_outcome,
    exposure = exposure,
    shade = FALSE,
    path = path_result_analysis,
    file_name = "plot_network_structure_outcome_exposure_non_shade_1"
)
