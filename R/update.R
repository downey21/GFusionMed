#' Internal function
#'
#' This is an internal function that is not exported.
updateUndirected <- function(v, dat,lambda,delta,Alpha,eta,kappa,pmat,no.de,bCb,no.tau) {
    # sampling alpha and kappa given eta
    # Input
    #   - v: current node
    #   - dat: nxp data matrix
    #   - Afull: chol2inv(chol((crossprod(dat,tdat)+diag(p)*(1/lambda))))
    #   - lambda, delta: hyperparameters for kappa
    #   - eta: pxp eta matrix
    #   - kappa: px1 kappa vector (inverse covariance for all nodes)
    #   - pmat: pxp matrix for hyper parameters for eta
    #   - no.de : number of directed edges connected to each node (for updating kappa) (vector)
    #   - no.de : t(b) %*% inv(C) %*% b for updating kappa for each node (vector)
    #   - no.tau:  number of nodes in the current layer

    n = nrow(dat)
    p = ncol(dat)

    # Graph movement
    if (sum(eta)==0) {
        is.swap = 0
    } else {
        is.swap = rbinom(1,1,1/2)
    }
    is.move = FALSE
    if (!is.swap == 1) { # ADD-DELETE
        w.ad = sample((1:p)[-v],1)
        new.eta = eta
        new.eta[v,w.ad] = new.eta[w.ad,v] = 1-eta[v,w.ad]
        w.up.l =  c(v,w.ad) # nodes for likelihood update
        is.move = TRUE
    } else { # SWAP
        cand0 = setdiff(which(eta[v,]==0),v)
        cand1 = which(eta[v,]==1)
        if (length(cand0)>0&length(cand1)>0) {
            w.ad1 = cand0[sample.int(length(cand0))[1]]
            w.ad2 = cand1[sample.int(length(cand1))[1]]
            new.eta = eta
            new.eta[v,w.ad1] = new.eta[w.ad1,v] = 1
            new.eta[v,w.ad2] = new.eta[w.ad2,v] = 0
            w.up.l = c(v,w.ad1,w.ad2)
            is.move = TRUE
        }
    }

    # Sampling parameters
    if (is.move) {

		Afull = base::chol2inv(base::chol((t_prod_cpp(dat,dat)+diag(p)*(1/lambda))))

		new.addr = sapply(w.up.l,function(x) which(new.eta[x,]==1),simplify=F)
		addr = sapply(w.up.l,function(x)which(eta[x,]==1),simplify=F)

		new.addr.inv = sapply(w.up.l,function(x) which(new.eta[x,]==0),simplify=F)
		addr.inv = sapply(w.up.l,function(x)which(eta[x,]==0),simplify=F)

		A.new = lapply(new.addr.inv,function(x) giveA(Afull,which=x))
		A = lapply(addr.inv,function(x) giveA(Afull,which=x))

		# eta
		lhr1 = sum(sapply(1:length(w.up.l),function(w) dmvnrm_arma(x = matrix(dat[,w.up.l[w]],ncol=n),mean=rep(0,n)
				,sigma = giveCov(dat[,new.addr[[w]],drop=F],rep(lambda,length(new.addr[[w]])),kappa[w.up.l[w]]),log=TRUE)))
		hyper1 = sum(sapply(w.up.l,function(w) sum(new.eta[w,-w]*log(pmat[w,-w]) + (1-new.eta[w,-w])*log(1-pmat[w,-w]))))
		lhr2 = sum(sapply(1:length(w.up.l),function(w) dmvnrm_arma(x = matrix(dat[,w.up.l[w]],ncol=n),mean=rep(0,n),
				sigma =giveCov(dat[,addr[[w]],drop=F],rep(lambda,length(addr[[w]])),kappa[w.up.l[w]]) ,log=TRUE)))
		hyper2 = sum(sapply(w.up.l,function(w) sum(eta[w,-w]*log(pmat[w,-w]) + (1-eta[w,-w])*log(1-pmat[w,-w]))))
		lhr = lhr1+hyper1-lhr2-hyper2

		if (log(runif(1))<lhr) {
			eta <- new.eta
			A <- A.new
		}

		for (i in 1:length(w.up.l)) {
			w = w.up.l[i]
			alpha = rep(0,p)
			resid = dat[,w]
			if (!is.null(A[[i]])){
			# alpha
			addr = which(eta[w,]==1)
				if (length(addr)>0) {
					alpha[addr] = rmvnrm_arma(1,prod_cpp(A[[i]],t_prod_cpp(dat[,addr,drop=F],dat[,w,drop=F])),A[[i]]/kappa[w])
				}
			Alpha[w,] = alpha
			# kappa
			resid = dat[,w] - prod_cpp(dat,as.matrix(alpha))
			}
			kappa[w] = stats::rgamma(1,shape=(n+delta+no.tau-1+no.de[w]+sum(alpha!=0))/2,rate=(lambda +sum(resid^2)+bCb[w]+(lambda +no.tau-1)*sum(alpha^2))/2)
		}
		
    }

    return(list(eta=eta,kappa=kappa,Alpha=Alpha,is.move=is.move))
}

#' Internal function
#'
#' This is an internal function that is not exported.
updateDirected <- function(y,X,lambda,delta,CCinv,B,Gamma,kappa,qmat,v.no.ue,v.aa,no.tau) {
	# sampling beta and kappa given gamma
	# Input
	#   - y: nx1 response vector corresponding to node v
	#   - X: nx(p*pP) covariate vector including parents of v and a function of parents of C(v)
	#   - lambda, delta: hyperparameters for kappa
	#   - CCinv: pxpP hyper parameter matrix for B
	#   - B : pxpP current B (permuted for response y)
	#   - Gamma: pxpP current gamma (permuted for response y)
	#   - kappa: a scalar for kappa for v (inverse covariance)
	#   - qmat: hyper parameter matrix for Gamma
	#   - v.no.ue: number of undirected edges for v (scalar) (for sampling kappa)
	#   - v.aa: t(alpha) %*% alpha (scalar) (for sampling kappa)
	#   - no.tau:  number of nodes in the current layer
    
	n = length(y)
    pP = ncol(B)
    p = nrow(B)
    vecpP = prod(dim(B))
    vqmat = c(qmat)

    # Graph movement
    if (sum(B[1,])==0) {
		is.swap=0
    } else {
		is.swap = rbinom(1,1,1/2)
	}
    is.move = FALSE
    if (!is.swap == 1) { # ADD-DELETE
		w.ad = sample(1:pP,1)
		new.Gamma = Gamma
		new.Gamma[1,w.ad] = 1-new.Gamma[1,w.ad]
		is.move = TRUE
    } else { # SWAP
		cand0 = which(Gamma[1,]==0)
		cand1 = which(Gamma[1,]==1)
		if (length(cand0)>0&length(cand1)>0) {
			w.ad1 = cand0[sample.int(length(cand0))[1]]
			w.ad2 = cand1[sample.int(length(cand1))[1]]
			new.Gamma = Gamma
			new.Gamma[1,w.ad1]  = 1
			new.Gamma[1,w.ad2]  = 0
			is.move = TRUE
		}
    }

    # Sampling parameters
    if (is.move) {
		vB = c(t(B))
		vGamma = c(t(Gamma))
		new.vGamma = c(t(new.Gamma))
		CC = 1/diag(CCinv)
		new.addr = which(new.vGamma==1)
		addr = which(vGamma==1)
		new.addr.inv = which(new.vGamma==0)
		addr.inv = which(vGamma==0)

		if (pP*p>n) {
			XC = X%*%diag(CC)
			tmp1 = base::chol2inv(base::chol(diag(n) + kappa * prod_t_cpp(XC,X)))
			Afull = diag(CC)-kappa*t_prod_cpp(XC,prod_cpp(tmp1,XC))
		} else {
			Afull = base::chol2inv(base::chol(kappa*t_prod_cpp(X,X)+CCinv))
		}

		A.new = giveA(Afull,which=new.addr.inv)
		A = giveA(Afull,which=addr.inv)

		# gamma
		lhr1 = dmvnrm_arma(x = matrix(y,ncol=n),mean=rep(0,n)
							   ,sigma = giveCov.dir(X=X[,new.addr,drop=F],CC=CC[new.addr],kappa=kappa),log=TRUE)
		hyper1 = sum(new.vGamma * log(vqmat) + (1-new.vGamma) * log(1-vqmat))
		lhr2 = dmvnrm_arma(x=matrix(y,ncol=n),mean=rep(0,n),sigma=giveCov.dir(X=X[,addr,drop=F],CC=CC[addr],kappa=kappa),log=TRUE)
		hyper2 = sum(vGamma * log(vqmat) + (1-vGamma) * log(1-vqmat))

		lhr = lhr1+hyper1-lhr2-hyper2

		if (log(runif(1))<lhr) {
			vGamma <- new.vGamma
			A <- A.new
		}

		# B
		vB = rep(0,length(vGamma))
		resid = y
		if (!is.null(A)) {
			ww = which(vGamma==1)
			if (length(ww)>0) {
				vB[ww]= rmvnrm_arma(1,kappa * prod_cpp(A,t_prod_cpp(X[,ww,drop=F],as.matrix(y))),A)
			}
			resid = y - prod_cpp(X,as.matrix(vB))
		}
		kappa = stats::rgamma(1,shape=(n+delta+no.tau-1+v.no.ue+sum(vB!=0))/2,rate=(lambda+sum(resid^2)+sum(vB^2/diag(CCinv))+(lambda +no.tau-1)*v.aa )/2)
		B = matrix(vB,ncol=pP,byrow=T)
		Gamma = matrix(vGamma,ncol=pP,byrow=T)
    }
	
    return(list(Gamma=Gamma,kappa=kappa,B=B,is.move=is.move))
}

#' Internal function
#'
#' This is an internal function that is not exported.
v.updateUndirected <- function(y,X,lambda,delta,Alpha,eta,kappa,pmat,no.de,bCb,no.tau) {
    # sampling alpha and kappa given eta
    # Input
    #   - y : response
    #   - X: design matrix
    #   - Afull: chol2inv(chol((crossprod(dat,tdat)+diag(p)*(1/lambda))))
    #   - lambda, delta: hyperparameters for kappa
    #   - eta: px1 eta matrix for v
    #   - kappa: px1 kappa vector (inverse covariance for all nodes)
    #   - pmat: pxp matrix for hyper parameters for eta
    #   - no.de : number of directed edges connected to each node (for updating kappa) (vector)
    #   - bCb : t(b) %*% inv(C) %*% b for updating kappa for each node (vector)
    #   - no.tau:  number of nodes in the current layer

    n = nrow(X)
    p = ncol(X)
    ### Graph movement ###
    if (sum(eta)==0) {is.swap=0
    }else{is.swap = rbinom(1,1,1/2)}
    is.move = FALSE
    if (!is.swap == 1) { #### ADD-DELETE
        w.ad = sample(1:p,1)
        new.eta = eta
        new.eta[w.ad] = 1-eta[w.ad]
        is.move = TRUE
    }else{ #### SWAP
        cand0 = which(eta==0)
        cand1 = which(eta==1)
        if (length(cand0)>0&length(cand1)>0) {
            w.ad1 = cand0[sample.int(length(cand0))[1]]
            w.ad2 = cand1[sample.int(length(cand1))[1]]
            new.eta = eta
            new.eta[w.ad1] = 1
            new.eta[w.ad2] = 0
            is.move = TRUE
        }
    }

    ### Sampling parameters ###
    if (is.move) {

        Afull = base::chol2inv(base::chol((crossprodCpp(X,X)+diag(p)*(1/lambda))))
        new.addr = which(new.eta==1)
        addr = which(eta==1)
        new.addr.inv =  which(new.eta==0)
        addr.inv = which(eta==0)
        A.new = giveA(Afull,which=new.addr.inv)
        A = giveA(Afull,which=addr.inv)

        #sample eta

        lhr1 = dmvnrm_arma(x = matrix(y,ncol=n),mean=rep(0,n),sigma = giveCov(X[,new.addr,drop=F],rep(lambda,length(new.addr)),kappa),log=TRUE)
        hyper1 = sum(new.eta*log(pmat) + (1-new.eta)*log(1-pmat))

        lhr2 = dmvnrm_arma(x = matrix(y,ncol=n),mean=rep(0,n),sigma =giveCov(X[,addr,drop=F],rep(lambda,length(addr)),kappa) ,log=TRUE)
        hyper2 = sum(eta*log(pmat) + (1-eta)*log(1-pmat))
        lhr = lhr1+hyper1-lhr2-hyper2


        if (log(runif(1))<lhr) {eta <- new.eta; A<-A.new}


        alpha = rep(0,p)
        resid = y
        if (!is.null(A)){
            # sample alpha
            addr = which(eta==1)
            if (length(addr)>0) {
                alpha[addr] = rmvnrm_arma(1,prodCpp(A,crossprodCpp(X[,addr,drop=F],as.matrix(y))),A/kappa)
            }
            # sample kappa
            resid = y -prodCpp(X,as.matrix(alpha))
        }
        kappa = stats::rgamma(1,shape=(n+delta+no.tau-1+no.de+sum(alpha!=0))/2,rate=(lambda +sum(resid^2)+bCb+(lambda +no.tau-1)*sum(alpha^2) )/2)

    } #if (is.move)
    return(list(eta=eta,kappa=kappa,Alpha=alpha,is.move=is.move))
}
