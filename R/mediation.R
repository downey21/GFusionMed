#' Perform Mediation Analysis
#'
#' This function prints.
#' @export
perform_mediation_analysis <- function(palist,chlist,fit.med,fit.out,ci=c(0.025,0.975)) {
    # - palist : parent list (used for fitting with the variable order)
    # - chlist : multilayer nodes
    # - fit.med: q-1 mediator model. (fit.med$Gamma includes no. mediator x (no. exposure + no. covariate)) (from second layer, does not allow GGM in the first layer)
    # - fit.out: outcome model
    # - ci: length 2 vector that includes upper and lower credible interval to calculate
    
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

    colnames(IEDout) = c("Exposure","Mediator","IED","CI_lower","CI_upper","PIP")
    colnames(IECout) = c("Exposure","Mediator","IEC","CI_lower","CI_upper","PIP")
    colnames(IEout) = c("Exposure","Mediator","IE","CI_lower","CI_upper","PIP")

    return(list(IED=IEDout,IEC=IECout,IE=IEout))
}