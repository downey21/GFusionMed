#' Mediation Analysis
#'
#' Conducts mediation analysis to identify indirect effects of omics variables on outcomes, 
#' revealing potential causal pathways between multi-omics data and the observed pharmacological response.
#'
#' @param result_structure The output from \code{fit_structure_model()} (Module 1 of GFusionMed).
#' @param result_outcome The output from \code{fit_outcome_model()} (Module 1 of GFusionMed).
#' @param exposure An optional character string specifying the exposure variable for 
#' which mediation effects are to be evaluated. If \code{NULL} (default), the function 
#' computes mediation effects for all possible exposure-mediator combinations.
#'
#' @details
#' This function performs mediation analysis using the outputs of \code{fit_structure_model()} 
#' and \code{fit_outcome_model()}. It computes indirect effects (IE) to uncover potential 
#' causal pathways between exposure and outcome through mediators. 
#' 
#' The analysis decomposes the total indirect effect (IE) into two components:
#'   \itemize{
#'     \item \code{IED}: The direct mediation effect between the exposure and mediator.
#'     \item \code{IEC}: The mediation effect arising from correlations within the same omics layer.
#'   }
#' Posterior inclusion probabilities (PIP) are computed using MCMC sampling to determine 
#' the significance of individual mediation effects. High PIP values for either IED or IEC 
#' suggest significant mediation effects.
#'
#' @return 
#' A data frame summarizing the mediation analysis results. The columns include:
#' \describe{
#'   \item{\code{Exposure}}{The exposure variable.}
#'   \item{\code{Mediator}}{The mediator variable.}
#'   \item{\code{Outcome}}{The outcome variable.}
#'   \item{\code{IED}}{The individual indirect effect between exposure and mediator.}
#'   \item{\code{PIP_IED}}{The posterior inclusion probability for IED.}
#'   \item{\code{IEC}}{The indirect effect from within-layer correlations.}
#'   \item{\code{PIP_IEC}}{The posterior inclusion probability for IEC.}
#'   \item{\code{IE}}{The individual indirect effect.}
#'   \item{\code{PIP_IE}}{The posterior inclusion probability for the total indirect effect.}
#' }
#'
#' @examples
#' # Mediation analysis
#' GFusionMed::perform_mediation_analysis(
#'   example_result_structure, example_result_outcome
#' )
#'
#' # Mediation analysis for a specific exposure variable
#' example_exposure <- "mRNA_EGFR"
#' 
#' GFusionMed::perform_mediation_analysis(
#'   example_result_structure, example_result_outcome, example_exposure
#' )
#' 
#' @references 
#' Seo, Dahun, et al. 
#' "Bayesian Multilayered Mediation Analysis for Cancer Pharmacogenomics." 
#' Stat 13.4 (2024)
#' 
#' @export
perform_mediation_analysis <- function(result_structure, result_outcome, exposure = NULL) {

    ci = c(0.025,0.975)

    if (length(result_structure) < 2) {
        stop("structure must contain at least 2 layers.")
    }

    fit.med <- list()
    fit.med <- result_structure[2:length(result_structure)]

    fit.out <- result_outcome

    information_outcome_model <- attr(result_outcome, "meta")

    palist <- information_outcome_model$palist
    chlist <- information_outcome_model$chlist
    nodenames <- information_outcome_model$column_name

    ov = do.call(c,chlist)
    layerinfo = do.call(c,sapply(1:length(chlist),function(k)rep(k,length(chlist[[k]]))))
    
    q = length(fit.med) + 1 # because the first layer model is not included
    p = length(ov) - 1
    
    inf.S = dim(result_structure[[1]]$eta)[3]
    
    if (q-1 == 1) {
        Aarray = c(sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
    } else {
        Aarray = do.call(c,sapply(1:(q-1),function(a) chlist[[a]],simplify=F))
    }

    v.m.getIndEffect <- function(v,m,palist,chlist,fit.med,fit.out,ci=c(0.025,0.975)){
        # - v : index for exposure
        # - m : position of layer of the whole mlGGM for fit.med
        # - palist : parent list (used for fitting with the variable order)
        # - chlist : multilayer nodes
        # - fit.med: mediator model. (fit.med$Gamma includes no. mediator x (no. exposure + no. covariate))
        # - fit.out: outcome model
        # - ci: length 2 vector that includes upper and lower credible interval to calculate
    
        v.med = chlist[[m]] # nodes in the mediator layer
        p.med = palist[[m]] # parents of the mediator layer m
        p.out = palist[[length(palist)]] # covariates (parents) of outcome model
        pos.v.med = which(p.med==v) ## position of v in parent set of mediator model (v is in the parent set)
        
        d.med = dim(fit.med$Gamma)
        d.out = dim(fit.out$Gamma)
        stopifnot(d.med[3]==d.out[2])
        
        inf.S = d.med[3]

        # IED
        outIED = outIEC = outIE = matrix(NA,nrow=length(v.med),ncol=4)
        outIE_inf_S = matrix(NA,nrow=length(v.med),ncol=inf.S)
        colnames(outIED) = c("IED","CI_lower","CI_upper","PIP")
        colnames(outIEC) = c("IEC","CI_lower","CI_upper","PIP")
        colnames(outIE) = c("IE","CI_lower","CI_upper","PIP")
        k = 0
        for (j in v.med) {
            k = k+1
            pos.j.med = which(v.med==j) # position of mediator j in mediator model
            pos.j.out = which(p.out==j) # position of mediator j in outcome model
            IED = sapply(1:inf.S,function(s) fit.out$B[pos.j.out,s] * fit.med$B[pos.j.med,pos.v.med,s] )
            IEC = sapply(1:inf.S,function(s) (-1) * fit.out$B[pos.j.out,s] * (fit.med$A[pos.j.med,-pos.j.med,s]%*%fit.med$B[-pos.j.med,pos.v.med,s]))
            IE = IED + IEC
            
            outIED[k,1] = mean(IED)
            outIED[k,2:3] = quantile(IED,probs=ci)
            outIED[k,4] = mean(IED!=0)
            
            outIEC[k,1] = mean(IEC)
            outIEC[k,2:3] = quantile(IEC,probs=ci)
            outIEC[k,4] = mean(IEC!=0)
            
            outIE[k,1] = mean(IE)
            outIE[k,2:3] = quantile(IE,probs=ci)
            outIE[k,4] = mean(IE!=0)

            outIE_inf_S[k,] = IE
        }

        return(list(IED=outIED,IEC=outIEC,IE=outIE,v.med=v.med,outIE_inf_S=outIE_inf_S))
    }

    if (!is.null(exposure)) {

        a <- which(nodenames[Aarray] == exposure)
        IEDout = IECout = IEout = numeric(0)

        a.layer= layerinfo[which(ov==a)]
        m.layer= (a.layer+1):q
        IED = IEC = IE = numeric(0)
        IE_inf_S = matrix(ncol = inf.S, nrow = 0)
        IE_layer = numeric(0)
        for (m in m.layer) {
            out = v.m.getIndEffect(v=a,m,palist,chlist,fit.med[[m-1]],fit.out,ci=c(0.025,0.975))
            IED = rbind(IED,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IED))
            IEC = rbind(IEC,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IEC))
            IE = rbind(IE,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IE))

            IE_layer_temp = apply(out$outIE_inf_S, 2, sum)
            IE_layer = rbind(IE_layer, c(m, mean(IE_layer_temp), mean(IE_layer_temp != 0)))

            IE_inf_S = rbind(IE_inf_S, out$outIE_inf_S)
        }

        rm(IE_layer_temp)

        IEDout = rbind(IEDout,IED)
        IECout = rbind(IECout,IEC)
        IEout = rbind(IEout,IE)

        IE_exposure_temp = apply(IE_inf_S, 2, sum)
        rm(IE_inf_S)
        IE_exposure = round(c(mean(IE_exposure_temp), mean(IE_exposure_temp != 0)), digits=8)

        DE_exposure_temp = result_outcome$B[a,1:inf.S]
        DE_exposure = round(c(mean(DE_exposure_temp), mean(DE_exposure_temp != 0)), digits=8 )

        TE_exposure_temp = DE_exposure_temp + IE_exposure_temp
        rm(IE_exposure_temp)
        rm(DE_exposure_temp)

        TE_exposure = round(c(mean(TE_exposure_temp), mean(TE_exposure_temp != 0)), digits=8)
        rm(TE_exposure_temp)

    } else {

        IEDout = IECout = IEout = numeric(0)
        for (a in Aarray) {
            a.layer= layerinfo[which(ov==a)]
            m.layer= (a.layer+1):q
            IED = IEC = IE = numeric(0)
            for (m in m.layer) {
                out = v.m.getIndEffect(v=a,m,palist,chlist,fit.med[[m-1]],fit.out,ci=c(0.025,0.975))
                IED = rbind(IED,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IED))
                IEC = rbind(IEC,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IEC))
                IE = rbind(IE,data.frame(rep(a,length(chlist[[m]])),out$v.med,out$IE))
            }
            IEDout = rbind(IEDout,IED)
            IECout = rbind(IECout,IEC)
            IEout = rbind(IEout,IE)
        }

    }

    colnames(IEDout) = c("Exposure","Mediator","IED","CI_lower","CI_upper","PIP")
    colnames(IECout) = c("Exposure","Mediator","IEC","CI_lower","CI_upper","PIP")
    colnames(IEout) = c("Exposure","Mediator","IE","CI_lower","CI_upper","PIP")

    IED <- as.matrix(IEDout)
    IEC <- as.matrix(IECout)
    IE <- as.matrix(IEout)

    w <- 1:nrow(IE)
    result <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(result) <- c("Exposure","Mediator","Outcome","individual_IE","PIP_individual_IE","IED","PIP_IED","IEC","PIP_IEC")
    if (length(w) > 0) {
        result <- matrix(NA, nrow = length(w), ncol = 9)
        result[ ,1:2] <- nodenames[IE[w, 1:2]]
        result[, 3] <- nodenames[length(nodenames)]
        result[ ,6] <- round(IED[w,"IED"], digit=8)
        result[ ,7] <- round(IED[w,"PIP"], digit=8)
        result[ ,8] <- round(IEC[w,"IEC"], digit=8)
        result[ ,9] <- round(IEC[w,"PIP"], digit=8)
        result[ ,4] <- round(IE[w,"IE"], digit=8)
        result[ ,5] <- round(IE[w,"PIP"], digit=8)

        colnames(result) <- c("Exposure","Mediator","Outcome","individual_IE","PIP_individual_IE","IED","PIP_IED","IEC","PIP_IEC")

        result <- as.data.frame(result)
    }

    if (!is.null(exposure)) {

        if (!(exposure %in% nodenames[1:(length(nodenames)-1)])) {
            stop("The exposure variable must be one of the variables in the structure.")
        }

        if (exposure %in% nodenames[chlist[[length(chlist) - 1]]]) {
            stop("Exposure variable cannot be in the layer immediately preceding the outcome layer.")
        }

        name_outcome <- nodenames[length(nodenames)]

        result_exposure <- result[result$Exposure == exposure, ]
        result_exposure <- result_exposure[, !(names(result_exposure) %in% c("Exposure", "Outcome"))]
        result_exposure <- result_exposure[order(result_exposure$PIP_individual_IE, decreasing = TRUE), ]

        rownames(result_exposure) = NULL

        colnames(IE_layer) = c("Layer", "layer_wise_IE", "PIP_layer_wise_IE")
        IE_layer = data.frame(IE_layer)

        IE_layer$Layer = attr(result_outcome, "meta")$structure_layer[IE_layer$Layer]

        attr(result_exposure, "info") <-
            list(
                Exposure = exposure,
                Outcome = name_outcome,
                TE = TE_exposure[1],
                PIP_TE = TE_exposure[2],
                DE = DE_exposure[1],
                PIP_DE = DE_exposure[2],
                IE = IE_exposure[1],
                PIP_IE = IE_exposure[2],
                layer_wise_IE = IE_layer,
                structure_layer = attr(result_outcome, "meta")$structure_layer,
                outcome_layer = attr(result_outcome, "meta")$outcome_layer
            )

        class(result_exposure) <- c("exposure_df", class(result_exposure))

        return(result_exposure)
    }

    return(result)
}

#' Internal function
#'
#' This is an internal function that is not exported.
.fmt_num <- function(x, digits = 6) {
    if (!is.numeric(x)) return(as.character(x))
    formatC(x, format = "f", digits = digits)
}

#' Internal function
#'
#' This is an internal function that is not exported.
.fmt_df_nums <- function(df, digits = 6) {
    if (!is.data.frame(df)) return(df)
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], .fmt_num, digits = digits)
    df
}

#' Internal function
#'
#' This is an internal function that is not exported.
.title_line <- function(x) cat(x, "\n", strrep("-", nchar(x)), "\n", sep = "")

#' Internal function
#'
#' This is an internal function that is not exported.
.kv <- function(label, value, digits = 6) {
    if (is.numeric(value)) value <- .fmt_num(value, digits)
    if (length(value) > 1) value <- paste(value, collapse = ", ")
    cat(sprintf("  %-14s: %s\n", label, value))
}

#' Internal function
#'
#' This is an internal function that is not exported.
.kv_pair <- function(label, val, pip, digits = 6) {
    if (is.numeric(val)) val <- .fmt_num(val, digits)
    if (is.numeric(pip)) pip <- .fmt_num(pip, digits)
    cat(sprintf("  %-14s: %s (%s)\n", paste0(label, " (PIP)"), val, pip))
}

#' Internal function
#'
#' This is an internal function that is not exported.
print.exposure_df <- function(x, ..., digits = 6, n_head = 10, show_all = FALSE) {
    info <- attr(x, "info")

    ## 1) Layer Info
    if (!is.null(info$structure_layer) || !is.null(info$outcome_layer)) {
        .title_line("Layer Info")
        if (!is.null(info$structure_layer)) .kv("structure_layer", info$structure_layer, digits)
        if (!is.null(info$outcome_layer))   .kv("outcome_layer",   info$outcome_layer, digits)
        cat("\n")
    }

    ## 2) Mediation Analysis
    if (!is.null(info)) {
        .title_line("Mediation Analysis")

        # Exposure / Outcome
        if (!is.null(info$Exposure)) .kv("Exposure", info$Exposure, digits)
        if (!is.null(info$Outcome))  .kv("Outcome",  info$Outcome,  digits)

        # TE/DE/IE (PIP)
        if (!is.null(info$TE) && !is.null(info$PIP_TE)) .kv_pair("TE", info$TE, info$PIP_TE, digits)
        if (!is.null(info$DE) && !is.null(info$PIP_DE)) .kv_pair("DE", info$DE, info$PIP_DE, digits)
        if (!is.null(info$IE) && !is.null(info$PIP_IE)) .kv_pair("IE", info$IE, info$PIP_IE, digits)

        cat("\n")
    }

    ## 3) Layer-wise IE Result
    if (!is.null(info$layer_wise_IE)) {
        .title_line("Layer-wise IE Result")
        lw <- info$layer_wise_IE
        lw <- .fmt_df_nums(as.data.frame(lw, stringsAsFactors = FALSE), digits = digits)

        print.data.frame(lw, row.names = FALSE)
        cat("\n")
    }

    ## 4) individual IE Result
    .title_line("individual IE Result")
    body_df <- if (isTRUE(show_all) || nrow(x) <= n_head) x else utils::head(x, n_head)

    body_df <- .fmt_df_nums(as.data.frame(body_df, stringsAsFactors = FALSE), digits = digits)
    print.data.frame(body_df, row.names = FALSE)
    if (!isTRUE(show_all) && nrow(x) > n_head) {
        cat(sprintf("... with %d more rows\n", nrow(x) - n_head))
    }

    invisible(x)
}
