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

#' Internal function
#'
#' This is an internal function that is not exported.
assign_colors <- function(vector_layer) {
  
    n <- length(vector_layer)
    palette_function <- colorRampPalette(RColorBrewer::brewer.pal(n = min(n, 8), name = "Set1"))
    colors <- palette_function(n)

    color_assignment <- setNames(colors, vector_layer)

    return(color_assignment)
}

#' Internal function
#'
#' This is an internal function that is not exported.
generate_x_centers <- function(vector_layer, total_range = 3000) {
    n <- length(vector_layer)
    max_position <- total_range / 2
    
    positions <- seq(from = -max_position, to = max_position, length.out = n)
    
    names(positions) <- vector_layer
    
    return(positions)
}