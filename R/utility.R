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

#' Internal function
#'
#' This is an internal function that is not exported.
assign_group <- function(node, layer_info) {
    
    node_split <- strsplit(node, "_")[[1]]
    layer <- node_split[1]
    node_name <- paste(node_split[-1], collapse = "_")
    
    if (layer %in% names(layer_info)) {
        
        group_info <- layer_info[[layer]]
        
        for (group_name in names(group_info)) {
            if (node_name %in% group_info[[group_name]]) {
                return(paste0(layer, "_", group_name))
            }
        }
    }
    
    stop("The group_information must contain the group information for every layer in the form of a list, ensuring no variables are missing from any layer.")
    
    return(NA)
}

#' Internal function
#'
#' This is an internal function that is not exported.
get_group_index <- function(group_name, node_index) {
    return(node_index[group_name])
}