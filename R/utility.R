#' Internal function
#'
#' This is an internal function that is not exported.
giveA <- function(Afull,which) { 
    if (length(which)>0) {
        rms::matinv(Afull,which=which,negate=F)[-which,-which,drop=F]
    } else {
        NULL
    }
}

#' Internal function
#'
#' This is an internal function that is not exported.
giveCov <- function(X,DD,kappa) {
	p = ncol(X)
	n = nrow(X)
	if (p>0) {
		(diag(n) + prod_t_cpp(prod_cpp(X,diag(DD,nrow=p)),X))/kappa
	} else {
        diag(n)/kappa
    }
}

#' Internal function
#'
#' This is an internal function that is not exported.
giveCov.dir <- function(X,CC,kappa) {
	p = ncol(X)
	n = nrow(X)
	if (p>0) {
		diag(n)/kappa + prod_t_cpp(prod_cpp(X,diag(CC,nrow=p)),X)
	} else {
        diag(n)/kappa
    }
}