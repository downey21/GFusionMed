#' Mediation Analysis
#'
#' Conducts mediation analysis to identify indirect effects of omics variables on outcomes, revealing potential causal pathways between multi-omics data and the observed pharmacological response.
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
    
    inf.S = ncol(fit.out$B)
    
    if (q-1 == 1) {
        A.array = c(sapply(1:(q-1),function(a) chlist[[a]],simplify=T))
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
        
        stopifnot()
        d.med = dim(fit.med$Gamma)
        d.out = dim(fit.out$Gamma)
        stopifnot(d.med[3]==d.out[2])
        
        inf.S = d.med[3]

        # IED
        outIED = outIEC = outIE = matrix(NA,nrow=length(v.med),ncol=4)
        colnames(outIED) = c("IED","CI_lower","CI_upper","PIP")
        colnames(outIEC) = c("IEC","CI_lower","CI_upper","PIP")
        colnames(outIE) = c("IE","CI_lower","CI_upper","PIP")
        k = 0
        for (j in v.med) {
            k = k+1
            pos.j.med = which(v.med==j) # position of mediator j in mediator model
            pos.j.out = which(p.out==j) # position of mediator j in outcome model
            IED = sapply(1:inf.S,function(s) fit.out$B[pos.j.out,s] * fit.med$B[pos.j.med,pos.v.med,s] )
            IEC = sapply(1:inf.S,function(s) fit.out$B[pos.j.out,s] * (fit.med$A[pos.j.med,-pos.j.med,s]%*%fit.med$B[-pos.j.med,pos.v.med,s]))
            IE = IED - IEC
            
            outIED[k,1] = mean(IED)
            outIED[k,2:3] = quantile(IED,probs=ci)
            outIED[k,4] = mean(IED!=0)
            
            outIEC[k,1] = mean(IEC)
            outIEC[k,2:3] = quantile(IEC,probs=ci)
            outIEC[k,4] = mean(IEC!=0)
            
            outIE[k,1] = mean(IE)
            outIE[k,2:3] = quantile(IE,probs=ci)
            outIE[k,4] = mean(IE!=0)
        }

        return(list(IED=outIED,IEC=outIEC,IE=outIE,v.med=v.med))
    }

    if (!is.null(exposure)) {

        a <- which(nodenames[Aarray] == exposure)
        IEDout = IECout = IEout = numeric(0)

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
    colnames(result) <- c("Exposure","Mediator","Outcome","IED","PIP_IED","IEC","PIP_IEC","IE","PIP_IE")
    if (length(w) > 0) {
        result <- matrix(NA, nrow = length(w), ncol = 9)
        result[ ,1:2] <- nodenames[IE[w, 1:2]]
        result[, 3] <- nodenames[length(nodenames)]
        result[ ,4] <- round(IED[w,"IED"], digit=8)
        result[ ,5] <- round(IED[w,"PIP"], digit=8)
        result[ ,6] <- round(IEC[w,"IEC"], digit=8)
        result[ ,7] <- round(IEC[w,"PIP"], digit=8)
        result[ ,8] <- round(IE[w,"IE"], digit=8)
        result[ ,9] <- round(IE[w,"PIP"], digit=8)

        colnames(result) <- c("Exposure","Mediator","Outcome","IED","PIP_IED","IEC","PIP_IEC","IE","PIP_IE")

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
        result_exposure <- result_exposure[order(result_exposure$PIP_IE, decreasing = TRUE), ]

        index_temp <- which(nodenames == exposure)

        DE <- round(rowMeans(result_outcome$B)[index_temp], digits = 8)
        PIP_DE <- round(rowMeans(result_outcome$Gamma)[index_temp], digits = 8)

        attr(result_exposure, "info") <-
            list(
                Exposure = exposure,
                Outcome = name_outcome,
                DE = DE,
                PIP_DE = PIP_DE,
                structure_layer = attr(result_outcome, "meta")$structure_layer,
                outcome_layer = attr(result_outcome, "meta")$outcome_layer
            )

        return(result_exposure)
    }

    return(result)
}


