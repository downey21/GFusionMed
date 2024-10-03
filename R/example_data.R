#' Example data for structure learning
#' 
#' Multi-omics (CNA, mRNA, Protein) data of lung cancer cell lines within the RTK pathway.
#' 
#' @format A list with 3 datasets
#' @source DepMap
"example_data_for_structure"

#' Example data for outcome learning
#' 
#' Multi-omics (CNA, mRNA, Protein) data of lung cancer cell lines within the RTK pathway, along with Erlotinib drug response (Drug) data.
#' 
#' @format A list with 4 datasets
#' @source DepMap
"example_data_for_outcome"

#' Example fitted result from structure learning
#' 
#' Output of GFusionMed::fit_structure_model(example_data_for_structure).
#' 
#' @format A list containing fitting information
"example_result_structure"

#' Example fitted result from outcome model learning
#' 
#' Output of GFusionMed::fit_outcome_model(data_for_outcome)).
#' 
#' @format A list containing fitting information
"example_result_outcome"
