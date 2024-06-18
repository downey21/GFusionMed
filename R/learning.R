#' Structure Learning
#'
#' This function prints the nu.
#' @export
fit_structure_model <- function(v.ch,v.pa,Y,eta.prob=0.1,gamma.prob=0.1,lambda,delta,burnin.S,inf.S) {
    # This function fits regression model for a chain component
    # Input
    # - v.ch : indices for the target chain component
    # - v.pa : indices for parents set
    # - Y : nxp data matrix
    stopifnot(length(v.ch)>2)
    S = burnin.S + inf.S
    dat.C = scale(Y[,v.ch,drop=F])
    n = nrow(dat.C)
    p = ncol(dat.C)
    pmat = matrix(0.2,ncol(dat.C),ncol(dat.C))
    diag(pmat) = 0

    if (is.null(v.pa)) {
        Gamma.list = NULL
        B.list = NULL
    } else {
        dat.P = scale(Y[,v.pa,drop=F])
        pP = ncol(dat.P)
        Gamma.list = B.list = array(0,dim=c(p,pP,inf.S))
        CC = matrix(1/lambda,ncol=pP,nrow=p) # hyper parameter for b for nonzero gamma
        qmat = matrix(0.1,ncol(dat.C),ncol(dat.P))
    }
    eta.list = A.list = array(0,dim=c(p,p,inf.S))
    kappa.list = matrix(0,nrow=inf.S,ncol=p)

    # initialization

    # eta (indicators for alpha)
    w.upper = which(upper.tri(diag(p)))
    eta = matrix(0,p,p)
    eta[w.upper] = rbinom(length(w.upper),size=1,prob=eta.prob)
    eta = eta + t_cpp(eta)
    diag(eta) = 0
    # Gamma (indicators for b)
    if (!is.null(v.pa)) {Gamma = matrix(rbinom(p*pP,size=1,prob=gamma.prob),nrow=p,ncol=pP)}
    # A (pxp alpha)
    A = pmat*eta
    # B (pxpP b)
    if (!is.null(v.pa)) {B = qmat*Gamma}
    # kappa (px1 vector)
    kappa = rep(1,p)

    s = 0

    while (s<S) {
        s = s+1
        if (s%%100==0) cat("Progress:",round(s/S*100, 1),"%","\n")
        for (v in sample(1:p)) { # in random order
            # Update eta, A, kappa
            if (!is.null(v.pa)) {
                tempDat = dat.C - prod_t_cpp(dat.P,B)
                no.de = rowSums(B!=0)
                bCb = sapply(1:p,function(v)sum(B[v,]^2/CC[v,]))
            } else {tempDat=dat.C
                no.de = bCb = rep(0,p)
            }
            up = updateUndirected(v=v,dat=tempDat,lambda=lambda,delta=delta,Alpha=A,eta=eta,kappa=kappa,pmat=pmat,no.de=no.de,bCb=bCb,no.tau=p)
            if (up$is.move) {
                eta = up$eta
                kappa = up$kappa
                A = up$Alpha
            }
            if (!is.null(v.pa)) {
                # Update Gamma, B, kappa
                v.ne = which(A[v,]!=0)
                v.ne.l = length(v.ne)
                vv.ne = c(v,v.ne)
                l.cl = length(vv.ne)
                if (v.ne.l>0) {
                    alpha = A[v,v.ne,drop=F]
                    y = c(dat.C[,v] - prod(dat.C[,v.ne,drop=F],alpha))
                    X = sapply(alpha,function(k) -k*dat.P,simplify=F)
                    X = cbind(dat.P,do.call(cbind,X))
                } else {
                    y = dat.C[,v]
                    X = dat.P
                    alpha = 0
                }
                CCinv = (1/CC[vv.ne,,drop=F]) * kappa[vv.ne]
                
                if (!all(dim(CCinv)==c(1,1))) {
                    if (ncol(CCinv)==1) {
                        CCinv = diag(c(CCinv))
                    } else {
                        CCinv = as.matrix(Matrix::bdiag(sapply(1:l.cl,function(x)diag(CCinv[x,]),simplify=F)))
                    }
                }
                
                tempB = B[vv.ne,,drop=F]
                tempGamma = Gamma[vv.ne,,drop=F]
                tempqmat = qmat[vv.ne,,drop=F]
                up = updateDirected(y=y,X=X,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa[v],qmat=tempqmat,v.no.ue =v.ne.l,v.aa=sum(alpha^2),no.tau=p)
                if (up$is.move){
                    Gamma[vv.ne,] = up$Gamma
                    kappa[v] = up$kappa
                    B[vv.ne,] = up$B
                }
            }
        }
        
        if (s>burnin.S) {
            ss = s - burnin.S
            if (!is.null(v.pa)) {
                Gamma.list[,,ss] = Gamma
                B.list[,,ss] = B
            }
            eta.list[,,ss] = eta
            A.list[,,ss] = A
            kappa.list[ss,] = kappa
        }
    }

    return(list(Gamma=Gamma.list,eta=eta.list,A=A.list,B=B.list,kappa=kappa.list))
}

#' Outcome Model Learning
#'
#' This function prints the nudd.
#' @export
fit_outcome_model <- function(y,X,gamma.prob=0.1,lambda,delta,burnin.S,inf.S,model="normal") {
    # Modified from ch.chaingraph function in the BANS package to allow one outcome
    # This function fits regression model for a chain component
    # Input
    # - y :nx1 vector for response (continuous or ordinal)
    # - X : nx(p+1) covariate matrix for exposure and mediators (ncol should be >2 with a exposure and a mediator)
    # - model: select "normal" or "probit"
    stopifnot(is.matrix(X),ncol(X)>=2)
    stopifnot(model %in% c("normal","probit"))
    n = nrow(X)
    S = burnin.S + inf.S
    
    if (model=="normal") {
        dat.C = scale(y)
    }else {
        zlevels = unique(y)
        K.z.levels = length(zlevels)
        z = as.numeric(factor(y,labels=1:K.z.levels))   ### Replace y to z
        gf = c(-Inf,stats::qnorm(1:(K.z.levels-1)/K.z.levels),Inf) ### g_0,....g_K Initialize g parameters
    }
    pmat = 0.2
    
    dat.P = scale(X)
    pP = ncol(dat.P)
    Gamma.list = B.list = matrix(NA,nrow=pP,ncol=inf.S)
    CC =rep(1/lambda,pP)# hyper parameter for b for nonzero gamma
    qmat = matrix(rep(0.1,ncol(dat.P)),ncol=pP)
    
    kappa.list = rep(0,nrow=inf.S)

    # initialization

    # Gamma (indicators for b)
    Gamma = matrix(stats::rbinom(pP,size=1,prob=gamma.prob),ncol=pP)
    # B (pxpP b)
    B = matrix(qmat*Gamma,ncol=pP)
    # kappa (px1 vector)
    kappa = 1
    
    s = 0

    while (s<S) {
        s = s+1
        if (s%%100==0) cat("Progress:",round(s/S*100, 1),"%","\n")
        CCinv = diag((1/CC) * kappa)
        tempB = B
        tempGamma = Gamma
        tempqmat = qmat
        
        if (model=="probit") {
            # Sample Y|Z,B,kappa  (Hoff and Niu p.212)
            ey = dat.P %*% c(tempB) # nx1 expected response for the latent
            a = gf[z] # nx1 lower bound
            b = gf[z+1] # n
            u = sapply(1:n,function(i) max(min(runif(1,pnorm((a[i]-ey[i])/kappa),pnorm((b[i]-ey[i])/kappa)),0.999),0.0001))
            y = ey + kappa * stats::qnorm(u) # Sampled y

            # Sample g_k|y,z, other g  (Hoff and Niu p.213) and Albert and Chib, 1993
            for (k in 2:K.z.levels) {
                a = max(max(y[z==(k-1)]),gf[k-1])
                b = max(min(min(y[z==k]),gf[k+1]),a+0.0001)
                gk = stats::runif(1,a,b)
                gf[k] = gk
            }            
            dat.C = y
        }
        
        up = updateDirected(y=dat.C,X=dat.P,lambda=lambda,delta=delta,CCinv=CCinv,B=tempB,Gamma=tempGamma,kappa=kappa,qmat=tempqmat,v.no.ue =0,v.aa=0,no.tau=1)
        if (up$is.move){
            Gamma[1,] = up$Gamma
            kappa = up$kappa
            B[1,] = up$B
        }

        if (s>burnin.S) {
            ss = s - burnin.S
            Gamma.list[,ss] = Gamma; B.list[,ss] = B
            kappa.list[ss] = kappa
        }
    }

    return(list(Gamma=Gamma.list,B=B.list,kappa=kappa.list))
}