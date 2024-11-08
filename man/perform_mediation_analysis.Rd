% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediation.R
\name{perform_mediation_analysis}
\alias{perform_mediation_analysis}
\title{Mediation Analysis}
\usage{
perform_mediation_analysis(result_structure, result_outcome, exposure = NULL)
}
\arguments{
\item{result_structure}{The output from \code{fit_structure_model()} (Module 1 of GFusionMed).}

\item{result_outcome}{The output from \code{fit_outcome_model()} (Module 1 of GFusionMed).}

\item{exposure}{An optional character string specifying the exposure variable for 
which mediation effects are to be evaluated. If \code{NULL} (default), the function 
computes mediation effects for all possible exposure-mediator combinations.}
}
\value{
A data frame summarizing the mediation analysis results. The columns include:
\describe{
  \item{\code{Exposure}}{The exposure variable.}
  \item{\code{Mediator}}{The mediator variable.}
  \item{\code{Outcome}}{The outcome variable.}
  \item{\code{IED}}{The individual indirect effect between exposure and mediator.}
  \item{\code{PIP_IED}}{The posterior inclusion probability for IED.}
  \item{\code{IEC}}{The indirect effect from within-layer correlations.}
  \item{\code{PIP_IEC}}{The posterior inclusion probability for IEC.}
  \item{\code{IE}}{The individual indirect effect.}
  \item{\code{PIP_IE}}{The posterior inclusion probability for the total indirect effect.}
}
}
\description{
Conducts mediation analysis to identify indirect effects of omics variables on outcomes, 
revealing potential causal pathways between multi-omics data and the observed pharmacological response.
}
\details{
This function performs mediation analysis using the outputs of \code{fit_structure_model()} 
and \code{fit_outcome_model()}. It computes indirect effects (IE) to uncover potential 
causal pathways between exposure and outcome through mediators. 

The analysis decomposes the total indirect effect (IE) into two components:
  \itemize{
    \item \code{IED}: The direct mediation effect between the exposure and mediator.
    \item \code{IEC}: The mediation effect arising from correlations within the same omics layer.
  }
Posterior inclusion probabilities (PIP) are computed using MCMC sampling to determine 
the significance of individual mediation effects. High PIP values for either IED or IEC 
suggest significant mediation effects.
}
\examples{
# Mediation analysis
GFusionMed::perform_mediation_analysis(
  example_result_structure, example_result_outcome
)

# Mediation analysis for a specific exposure variable
example_exposure <- "mRNA_EGFR"

GFusionMed::perform_mediation_analysis(
  example_result_structure, example_result_outcome, example_exposure
)

}
\references{
Seo, Dahun, et al. 
"Bayesian Multilayered Mediation Analysis for Cancer Pharmacogenomics." 
Stat 13.4 (2024)
}
