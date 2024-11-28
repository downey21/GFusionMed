#' Summarize the edge information from fitted results
#'
#' This function calculates the edge strength from the MCMC sampling results, 
#' providing a summary of both undirected and directed edges.
#'
#' @param x The output from either \code{fit_structure_model()} or 
#' \code{fit_outcome_model()} from Module 1 of the GFusionMed package.
#'
#' @details 
#' The function summarizes both undirected and directed edges based on the MCMC 
#' sampling results. It returns a data frame with the start and end nodes of each edge, 
#' the edge type (either "und" for undirected or "dir" for directed), and the edge strength, 
#' which reflects the frequency or strength of the inferred relationship across MCMC samples.
#'
#' @return 
#' A data frame with the following columns:
#' \describe{
#'   \item{\code{Start}}{The starting node of the edge.}
#'   \item{\code{End}}{The ending node of the edge.}
#'   \item{\code{Type}}{The type of edge, either "und" (undirected) or "dir" (directed).}
#'   \item{\code{Strength}}{The strength of the edge, calculated based on MCMC samples.}
#' }
#'
#' @examples
#' data(example_result_structure, package = "GFusionMed")
#' data(example_result_outcome, package = "GFusionMed")
#'
#' GFusionMed::edge_summary(example_result_structure)
#' GFusionMed::edge_summary(example_result_outcome)
#' 
#' @export
edge_summary <- function(x, ...) {
    UseMethod("edge_summary")
}

#' Summarize the edge information from fitted result
#'
#' This function calculates edge strength from the MCMC sampling results.
#' @export
edge_summary.structure <- function(result_structure, exposure = NULL) {
    information_structure_model <- attr(result_structure, "meta")

    palist <- information_structure_model$palist
    chlist <- information_structure_model$chlist
    nodenames <- information_structure_model$column_name
    vector_layer <- attr(result_structure, "meta")$structure_layer

    p <- length(nodenames)
    q <- length(chlist)
    
    B = matrix(0,p,p)
    cB = numeric(0)
    for (j in 1:q) {
        tmp = apply(result_structure[[j]]$eta,c(1,2),mean)
        B[chlist[[j]],chlist[[j]]] = tmp
        cB = c(cB,tmp[lower.tri(tmp)])
        if (j!=1) {
            tmp = c(apply(result_structure[[j]]$Gamma,c(1,2),mean))
            B[as.matrix(expand.grid(chlist[[j]],palist[[j]]))] = tmp
            cB = c(cB,tmp)
        }
    }
    o = order(cB,decreasing=T)
    index = which(cumsum(1-cB[o]) * (1/1:length(cB)) < 0.1) # FDR at 0.1
    ii = ifelse(length(index) == 0, NA, max(index))
    cut = cB[o][ii]
    
    matB = matrix(0,p,p)

    # matB[B>=cut] = 1
    matB[B>0] = 1
    
    und.w = which(matB + t(matB)==2,arr.ind=T)
    dir.w = which(matB==1 & t(matB)==0,arr.ind=T)
    
    if (length(und.w) > 0) {
        und.w = unique(t(apply(und.w,1,sort)))
        undattr = cbind(nodenames[und.w[,1]],nodenames[und.w[,2]],"und",round(B[und.w],digit=8))
    } else {
        undattr = NULL
    }
    
    if (length(dir.w) > 0) {
        dirattr = cbind(nodenames[dir.w[,2]],nodenames[dir.w[,1]],"dir",round(B[dir.w],digit=8))
    } else {
        dirattr = NULL
    }
    
    edgeattr = rbind(undattr,dirattr)

    colnames(edgeattr) <- c("Start", "End", "Type", "Strength")
    edgeattr <- as.data.frame(edgeattr)

    edgeattr$Strength <- as.numeric(edgeattr$Strength)

    if (!is.null(exposure)) {

        if (!(exposure %in% nodenames)) {
            stop("The exposure variable must be one of the variables in the structure.")
        }

        if (exposure %in% nodenames[chlist[[length(chlist)]]]) {
            stop("The exposure variable cannot be in the last layer of the structure.")
        }

        exposure_layer_part <- strsplit(exposure, "_")[[1]][1]
        exposure_layer <- which(vector_layer == exposure_layer_part)

        if (exposure_layer >= 2) {

            patterns_to_remove <- paste0("^", vector_layer[1:(exposure_layer - 1)], "_")

            for (pattern in patterns_to_remove) {
                edgeattr <- edgeattr[!grepl(pattern, edgeattr$Start), ]
            }
        
        }

        pattern_for_layer <- paste0("^", vector_layer[exposure_layer], "_")
    
        edgeattr <- edgeattr[!(grepl(pattern_for_layer, edgeattr$Start) & edgeattr$Start != exposure), ]

        edgeattr <- edgeattr[!(grepl(pattern_for_layer, edgeattr$End)), ]

        attr(edgeattr, "info") <-
            list(
                Exposure = exposure,
                structure_layer = vector_layer
            )

    }

    if (is.null(exposure)) {

        attr(edgeattr, "info") <-
            list(
                structure_layer = vector_layer
            )

    }

    return(edgeattr)
}

#' Summarize the edge information from fitted result
#'
#' This function calculates edge strength from the MCMC sampling results.
#' @export
edge_summary.outcome <- function(result_outcome, exposure = NULL) {
    information_outcome_model <- attr(result_outcome, "meta")

    palist <- information_outcome_model$palist
    chlist <- information_outcome_model$chlist
    nodenames <- information_outcome_model$column_name

    cB = rowMeans(result_outcome$Gamma)
    o = order(cB,decreasing=T)
    index = which(cumsum(1-cB[o]) * (1/1:length(cB)) < 0.1) # FDR at 0.1
    ii = ifelse(length(index) == 0, NA, max(index))
    cut = cB[o][ii]

    # w = which(cB>=cut)
    # w = 1:length(cB)
    
    if (!is.null(exposure)) {

        if (!(exposure %in% nodenames[1:(length(nodenames)-1)])) {
            stop("The exposure variable must be one of the variables in the structure.")
        }

        if (exposure %in% nodenames[chlist[[length(chlist) - 1]]]) {
            stop("Exposure variable cannot be in the layer immediately preceding the outcome layer.")
        }

        w <- which(nodenames == exposure)

        dname = nodenames[chlist[[4]]]
        ddir = cbind(nodenames[w],rep(dname,length(w)),rep("dir",length(w)),cB[w])

        colnames(ddir) <- c("Start", "End", "Type", "Strength")
        ddir <- as.data.frame(ddir)

        ddir$Strength <- as.numeric(ddir$Strength)

        attr(ddir, "info") <-
            list(
                Exposure = exposure,
                Outcome = nodenames[length(nodenames)],
                structure_layer = attr(result_outcome, "meta")$structure_layer,
                outcome_layer = attr(result_outcome, "meta")$outcome_layer
            )

    } else {

        w = 1:length(cB)
        
        if (length(w)>0) {
            dname = nodenames[chlist[[4]]]
            ddir = cbind(nodenames[w],rep(dname,length(w)),rep("dir",length(w)),cB[w])
        }
        
        colnames(ddir) <- c("Start", "End", "Type", "Strength")
        ddir <- as.data.frame(ddir)

        ddir$Strength <- as.numeric(ddir$Strength)

        attr(ddir, "info") <-
            list(
                structure_layer = attr(result_outcome, "meta")$structure_layer,
                outcome_layer = attr(result_outcome, "meta")$outcome_layer
            )

    }

    return(ddir)
}

#' Network Visualization
#'
#' Visualizes the network, where edge transparency reflects the strength of inferred relationships.
#' Depending on the input, the visualization can focus on either the multi-layered structure, 
#' the outcome model, or both. It also supports focused visualization for a specific exposure variable.
#'
#' @param result_structure The output from \code{fit_structure_model()} (Module 1 of GFusionMed). 
#' If provided, the network visualization will reflect the learned multi-layered structure.
#' @param result_outcome The output from \code{fit_outcome_model()} (Module 1 of GFusionMed). 
#' If provided, the network visualization will reflect the learned outcome model.
#' @param exposure An optional character string specifying the exposure variable to focus on. 
#' If \code{NULL} (default), the entire network is visualized. When a specific exposure is provided, 
#' the visualization will focus on the pathway from the exposure to the outcome through mediators.
#' @param shade A logical value indicating whether to dim unrelated variables when focusing on 
#' a specific exposure. Defaults to \code{TRUE}, meaning unrelated variables are dimmed. If \code{FALSE}, 
#' unrelated variables are excluded from the visualization entirely. This parameter only applies when 
#' \code{exposure} is not \code{NULL}.
#' @param path A character string specifying the directory path where the output PDF file will be saved.
#' @param file_name A character string specifying the name of the output PDF file (without the ".pdf").
#' @param width A numeric value specifying the width of the PDF file. Default is 20.
#' @param height A numeric value specifying the height of the PDF file. Default is 11.
#' @param vector_layer_color An optional vector specifying the colors for each layer in the network. 
#' If not provided, default colors will be used.
#' @param seed An integer value for setting the random seed to ensure reproducibility in node placement.
#'
#' @details 
#' The \code{plot_network()} function creates PDF visualizations of the network, where the edge transparency 
#' indicates the strength of the inferred relationships. It supports visualizing:
#' \itemize{
#'   \item The multi-layered structure network when only \code{result_structure} is provided.
#'   \item The outcome model network when only \code{result_outcome} is provided.
#'   \item A combined network of the structure and outcome model when both inputs are provided.
#' }
#' If a specific exposure variable is provided, the visualization focuses on the pathway from the exposure 
#' to the outcome through mediators. The \code{shade} parameter controls whether unrelated variables are 
#' dimmed (\code{TRUE}) or excluded (\code{FALSE}).
#'
#' The output is saved as a PDF file at the specified \code{path} and with the specified \code{file_name}.
#' No value is returned by the function.
#'
#' @return 
#' This function does not return a value. The visualization is saved as a PDF file.
#'
#' @examples
#' # Structure network visualization
#' GFusionMed::plot_network(
#'   result_structure = example_result_structure,
#'   path = "Set your path", file_name = "plot_network_structure"
#' )
#'
#' # Outcome network visualization
#' GFusionMed::plot_network(
#'   result_outcome = example_result_outcome,
#'   path = "Set your path", file_name = "plot_network_outcome"
#' )
#'
#' # Structure-outcome network visualization
#' GFusionMed::plot_network(
#'   result_structure = example_result_structure,
#'   result_outcome = example_result_outcome, 
#'   path = "Set your path", file_name = "plot_network_structure_outcome"
#' )
#'
#' # Network visualization for a specific exposure variable
#' example_exposure <- "mRNA_EGFR"
#' GFusionMed::plot_network(
#'   result_structure = example_result_structure,
#'   result_outcome = example_result_outcome, 
#'   exposure = example_exposure,
#'   path = "Set your path", file_name = "plot_network_structure_outcome_exposure"
#' )
#'
#' @export
plot_network <- function(result_structure = NULL, result_outcome = NULL, exposure = NULL, shade = TRUE, path, file_name, width = 20, height = 11, vector_layer_color = NULL, seed = 1234) {

    if (is.null(result_structure) & !is.null(result_outcome)) {
        return(
            plot_network.outcome(
                result_outcome = result_outcome,
                path = path,
                file_name = file_name,
                width = width,
                height = height,
                vector_layer_color = vector_layer_color,
                seed = seed
            )
        )
    } else if (!is.null(result_structure) & is.null(result_outcome)) {
        return(
            plot_network.structure(
                result_structure = result_structure,
                path = path,
                file_name = file_name,
                width = width,
                height = height,
                vector_layer_color = vector_layer_color,
                seed = seed
            )
        )
    } else if (!is.null(result_structure) & !is.null(result_outcome)) {
        return(
            plot_network.structure_outcome(
                result_structure = result_structure,
                result_outcome = result_outcome,
                exposure = exposure,
                shade = shade,
                path = path,
                file_name = file_name,
                width = width,
                height = height,
                vector_layer_color = vector_layer_color,
                seed = seed
            )
        )
    } else if (is.null(result_structure) & is.null(result_outcome)) {
        stop("At least one of the two must not be NULL")
    }

}

#' Internal function
#'
#' This is an internal function that is not exported.
plot_network.structure <- function(result_structure, path, file_name, width = 20, height = 11, vector_layer_color = NULL, seed = 1234) {

    data_edge <- edge_summary(result_structure)

    set.seed(seed) 

    vector_layer <- attr(data_edge, "info")$structure_layer

    if (!is.null(vector_layer_color)) {

        if (!is.vector(vector_layer_color)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of structure layers.")
        }

        if (length(vector_layer_color) != length(vector_layer)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of structure layers.")
        }

    }

    if (!is.null(vector_layer_color)) {
        vector_layer_color <- setNames(vector_layer_color, vector_layer)
    }

    if (is.null(vector_layer_color)) {
        vector_layer_color <- assign_colors(vector_layer)
    }

    node_layer <-
        sapply(
            unique(c(data_edge$Start, data_edge$End)),
            function(x) strsplit(x, "_")[[1]][1]
        )

    node_colors <- setNames(vector_layer_color[node_layer], names(node_layer))

    g <- igraph::graph_from_data_frame(data_edge, directed = TRUE)
    igraph::V(g)$color <- node_colors[igraph::V(g)$name]
    igraph::V(g)$label <- sub("^[^_]+_", "", igraph::V(g)$name)
    igraph::E(g)$weight <- data_edge$Strength^2
    igraph::E(g)$lty <- 1
    igraph::E(g)$arrow.size <- ifelse(data_edge$Type == "dir", 2.5, 3.5)
    igraph::E(g)$arrow.mode <- ifelse(data_edge$Type == "dir", 2, 0)
    igraph::E(g)$color <- apply(as.matrix(igraph::E(g)$weight), 1, function(w) {
        rgb(0.5, 0.5, 0.5, alpha = w / max(igraph::E(g)$weight))
    })

    node_positions <- matrix(NA, nrow = length(igraph::V(g)), ncol = 2)

    R <- 2
    x_centers <- generate_x_centers(vector_layer)
    y_centers <- stats::runif(length(vector_layer), -120, 120)
    names(y_centers) <- names(x_centers)

    for (category in names(x_centers)) {
        group_nodes <- which(node_layer == category)
        num_nodes <- length(group_nodes)
        theta_step <- 2 * pi / num_nodes
        
        for (i in seq_along(group_nodes)) {
            r <- sqrt(runif(1, 15000, 15000.5)) * R
            theta <- theta_step * i
            
            node_positions[group_nodes[i], 1] <- x_centers[category] + r * cos(theta) + stats::runif(1, -120, 120)
            node_positions[group_nodes[i], 2] <- y_centers[category] + r * sin(theta) + stats::runif(1, -120, 120)
        }
    }

    pdf(file.path(path, paste0(file_name, ".pdf")), width = width, height = height)

    igraph::plot.igraph(g,
                        edge.width = 3,
                        vertex.size = 15,
                        vertex.label.cex = 1.6,
                        vertex.label.color = "black",
                        edge.arrow.size = igraph::E(g)$arrow.size,
                        edge.arrow.mode = igraph::E(g)$arrow.mode,
                        edge.color = igraph::E(g)$color,
                        layout = node_positions,
                        asp = 0
    )

    legend(
        "bottomright", 
        legend = names(vector_layer_color),
        col = vector_layer_color,
        pch = 21,
        pt.bg = vector_layer_color,
        pt.cex = 3,
        cex = 1.6,
        bty = "n"
    )

    invisible(dev.off())

}

#' Internal function
#'
#' This is an internal function that is not exported.
plot_network.outcome <- function(result_outcome, path, file_name, width = 20, height = 11, vector_layer_color = NULL, seed = 1234) {

    data_edge <- edge_summary(result_outcome)
    
    set.seed(seed) 

    vector_layer <- c(attr(data_edge, "info")$structure_layer, attr(data_edge, "info")$outcome_layer)

    if (!is.null(vector_layer_color)) {

        if (!is.vector(vector_layer_color)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of structure layers.")
        }

        if (length(vector_layer_color) != length(vector_layer)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of structure layers.")
        }

    }

    if (!is.null(vector_layer_color)) {
        vector_layer_color <- setNames(vector_layer_color, vector_layer)
    }

    if (is.null(vector_layer_color)) {
        vector_layer_color <- assign_colors(vector_layer)
    }

    node_layer <-
        sapply(
            unique(c(data_edge$Start, data_edge$End)),
            function(x) strsplit(x, "_")[[1]][1]
        )

    node_colors <- setNames(vector_layer_color[node_layer], names(node_layer))

    g <- igraph::graph_from_data_frame(data_edge, directed = TRUE)
    igraph::V(g)$color <- node_colors[igraph::V(g)$name]
    igraph::V(g)$label <- sub("^[^_]+_", "", igraph::V(g)$name)
    igraph::E(g)$weight <- data_edge$Strength^2
    igraph::E(g)$lty <- 1
    igraph::E(g)$arrow.size <- ifelse(data_edge$Type == "dir", 2.5, 3.5)
    igraph::E(g)$arrow.mode <- ifelse(data_edge$Type == "dir", 2, 0)
    igraph::E(g)$color <- apply(as.matrix(igraph::E(g)$weight), 1, function(w) {
        rgb(0.5, 0.5, 0.5, alpha = w / max(igraph::E(g)$weight))
    })

    node_positions <- matrix(NA, nrow = length(igraph::V(g)), ncol = 2)

    R <- 2
    x_centers <- generate_x_centers(vector_layer)
    y_centers <- stats::runif(length(vector_layer), -120, 120)
    names(y_centers) <- names(x_centers)

    for (category in names(x_centers)) {
        group_nodes <- which(node_layer == category)
        num_nodes <- length(group_nodes)
        theta_step <- 2 * pi / num_nodes
        
        for (i in seq_along(group_nodes)) {
            r <- sqrt(runif(1, 15000, 15000.5)) * R
            theta <- theta_step * i
            
            node_positions[group_nodes[i], 1] <- x_centers[category] + r * cos(theta) + stats::runif(1, -120, 120)
            node_positions[group_nodes[i], 2] <- y_centers[category] + r * sin(theta) + stats::runif(1, -120, 120)
        }
    }

    pdf(file.path(path, paste0(file_name, ".pdf")), width = width, height = height)

    igraph::plot.igraph(g,
                        edge.width = 3,
                        vertex.size = 15,
                        vertex.label.cex = 1.6,
                        vertex.label.color = "black",
                        edge.arrow.size = igraph::E(g)$arrow.size,
                        edge.arrow.mode = igraph::E(g)$arrow.mode,
                        edge.color = igraph::E(g)$color,
                        layout = node_positions,
                        asp = 0
    )

    legend(
        "bottomright", 
        legend = names(vector_layer_color),
        col = vector_layer_color,
        pch = 21,
        pt.bg = vector_layer_color,
        pt.cex = 3,
        cex = 1.6,
        bty = "n"
    )

    invisible(dev.off())

}

#' Internal function
#'
#' This is an internal function that is not exported.
plot_network.structure_outcome <- function(result_structure, result_outcome, exposure = NULL, shade = TRUE, path, file_name, width = 20, height = 11, vector_layer_color = NULL, seed = 1234) {

    set.seed(seed) 

    if (!is.null(exposure)) {

        if (shade) {

            data_edge <-
                rbind(
                    GFusionMed::edge_summary(result_outcome),
                    GFusionMed::edge_summary(result_structure)
                )

        } else {

            data_edge <-
                rbind(
                    GFusionMed::edge_summary(result_outcome),
                    GFusionMed::edge_summary(result_structure, exposure)
                )

        }

    } else {

        data_edge <-
            rbind(
                GFusionMed::edge_summary(result_outcome),
                GFusionMed::edge_summary(result_structure)
            )

    }

    vector_layer <- c(attr(data_edge, "info")$structure_layer, attr(data_edge, "info")$outcome_layer)

    if (!is.null(exposure) && shade == FALSE) {

        exposure_layer_part <- strsplit(exposure, "_")[[1]][1]
        exposure_layer <- which(vector_layer == exposure_layer_part)

        if (exposure_layer >= 2) {

            patterns_to_remove <- paste0("^", vector_layer[1:(exposure_layer - 1)], "_")

            for (pattern in patterns_to_remove) {
                data_edge <- data_edge[!grepl(pattern, data_edge$Start), ]
            }
        
        }

        pattern_for_layer <- paste0("^", vector_layer[exposure_layer], "_")
    
        data_edge <- data_edge[!(grepl(pattern_for_layer, data_edge$Start) & data_edge$Start != exposure), ]

    }

    if (!is.null(vector_layer_color)) {

        if (!is.vector(vector_layer_color)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of layers, including the outcome layer.")
        }

        if (length(vector_layer_color) != length(vector_layer)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of layers, including the outcome layer.")
        }

    }

    if (!is.null(vector_layer_color)) {
        vector_layer_color <- setNames(vector_layer_color, vector_layer)
    }

    if (is.null(vector_layer_color)) {
        vector_layer_color <- assign_colors(vector_layer)
    }

    if (!is.null(exposure) && shade == FALSE) {
        vector_layer_color_for_legend <- vector_layer_color[exposure_layer:length(vector_layer_color)]
    }

    node_layer <-
        sapply(
            unique(c(data_edge$Start, data_edge$End)),
            function(x) strsplit(x, "_")[[1]][1]
        )

    node_colors <- setNames(vector_layer_color[node_layer], names(node_layer))

    if (!is.null(exposure) && shade == TRUE) {

        exposure_layer_part <- strsplit(exposure, "_")[[1]][1]
        exposure_layer <- which(vector_layer == exposure_layer_part)

        if (exposure_layer >= 2) {

            patterns_to_remove <- paste0("^", vector_layer[1:(exposure_layer - 1)], "_")

            for (name in names(node_colors)) {
                for (pattern in patterns_to_remove) {
                    if (grepl(pattern, name)) {
                        node_colors[name] <- adjustcolor(node_colors[name], alpha.f = 0.25)
                    }
                }
            }
        
            for (pattern in patterns_to_remove) {
                data_edge[grepl(pattern, data_edge$Start), "Strength"] <- data_edge[grepl(pattern, data_edge$Start), "Strength"] / 4
            }

        }

        pattern_for_layer <- paste0("^", vector_layer[exposure_layer], "_")
    
        for (name in names(node_colors)) {
            if (grepl(pattern_for_layer, name) & name != exposure) {
                node_colors[name] <- adjustcolor(node_colors[name], alpha.f = 0.25)
            }
        }

        data_edge[grepl(pattern_for_layer, data_edge$Start) & data_edge$Start != exposure, "Strength"] <- data_edge[grepl(pattern_for_layer, data_edge$Start) & data_edge$Start != exposure, "Strength"] / 4
        data_edge[grepl(pattern_for_layer, data_edge$Start) & grepl(pattern_for_layer, data_edge$End), "Strength"] <- data_edge[grepl(pattern_for_layer, data_edge$Start) & grepl(pattern_for_layer, data_edge$End), "Strength"] / 4

    }

    g <- igraph::graph_from_data_frame(data_edge, directed = TRUE)
    igraph::V(g)$color <- node_colors[igraph::V(g)$name]
    igraph::V(g)$label <- sub("^[^_]+_", "", igraph::V(g)$name)
    igraph::E(g)$weight <- data_edge$Strength^2
    igraph::E(g)$lty <- 1
    igraph::E(g)$arrow.size <- ifelse(data_edge$Type == "dir", 2.5, 3.5)
    igraph::E(g)$arrow.mode <- ifelse(data_edge$Type == "dir", 2, 0)
    igraph::E(g)$color <- apply(as.matrix(igraph::E(g)$weight), 1, function(w) {
        rgb(0.5, 0.5, 0.5, alpha = w / max(igraph::E(g)$weight))
    })

    node_positions <- matrix(NA, nrow = length(igraph::V(g)), ncol = 2)

    R <- 2
    x_centers <- generate_x_centers(vector_layer)
    y_centers <- stats::runif(length(vector_layer), -120, 120)
    names(y_centers) <- names(x_centers)

    for (category in names(x_centers)) {
        group_nodes <- which(node_layer == category)
        num_nodes <- length(group_nodes)
        theta_step <- 2 * pi / num_nodes
        
        for (i in seq_along(group_nodes)) {
            r <- sqrt(runif(1, 15000, 15000.5)) * R
            theta <- theta_step * i
            
            node_positions[group_nodes[i], 1] <- x_centers[category] + r * cos(theta) + stats::runif(1, -120, 120)
            node_positions[group_nodes[i], 2] <- y_centers[category] + r * sin(theta) + stats::runif(1, -120, 120)
        }
    }

    pdf(file.path(path, paste0(file_name, ".pdf")), width = width, height = height)

    igraph::plot.igraph(g,
                        edge.width = 3,
                        vertex.size = 15,
                        vertex.label.cex = 1.6,
                        vertex.label.color = "black",
                        edge.arrow.size = igraph::E(g)$arrow.size,
                        edge.arrow.mode = igraph::E(g)$arrow.mode,
                        edge.color = igraph::E(g)$color,
                        layout = node_positions,
                        asp = 0
    )

    if (!is.null(exposure) && shade == FALSE) {
        legend(
            "bottomright", 
            legend = names(vector_layer_color_for_legend),
            col = vector_layer_color_for_legend,
            pch = 21,
            pt.bg = vector_layer_color_for_legend,
            pt.cex = 3,
            cex = 1.6,
            bty = "n"
        )    
    } else {
        legend(
            "bottomright", 
            legend = names(vector_layer_color),
            col = vector_layer_color,
            pch = 21,
            pt.bg = vector_layer_color,
            pt.cex = 3,
            cex = 1.6,
            bty = "n"
        )
    }

    invisible(dev.off())

}

#' Sankey Diagram
#'
#' Visualizes the network using a Sankey diagram when group information for the variables is available.
#' The diagram illustrates the flow and relationships between groups of variables across layers, 
#' providing a clear visualization of how variables and outcomes are connected.
#'
#' @param edge_structure The output from \code{GFusionMed::edge_summary()} using 
#' the results from \code{fit_structure_model()} (Module 1 of GFusionMed).
#' @param edge_outcome The output from \code{GFusionMed::edge_summary()} using 
#' the results from \code{fit_outcome_model()}. If multiple outcome results are available, 
#' they can be combined using \code{rbind()}.
#' @param group_information A list containing the group assignments for each layer. 
#' All variables must belong to a group, and no variable can be left unassigned.
#' @param vector_layer_color An optional vector specifying the colors for each layer in the network. 
#' If not provided, default colors will be used.
#'
#' @details 
#' The \code{plot_sankey()} function generates a Sankey diagram to visualize relationships across 
#' layers or groups of variables. It integrates the structure and outcome networks to display how variables interact across layers.
#' Group information is required for the variables in each layer, and any missing group assignment will result in an error.
#'
#' @return 
#' This function outputs a Sankey diagram to the plotting device.
#'
#' @examples
#' # Generate edge summaries
#' edge_structure <- GFusionMed::edge_summary(example_result_structure)
#' edge_outcome <- GFusionMed::edge_summary(example_result_outcome)
#'
#' # Define group information for variables
#' group_information <- list(
#'   CNA = list(
#'     group1 = c("SHC1", "ERBB3", "ERBB2"),
#'     group2 = c("SRC", "EGFR")
#'   ),
#'   mRNA = list(
#'     group1 = c("SHC1", "ERBB3", "ERBB2"),
#'     group2 = c("SRC", "EGFR")
#'   ),
#'   Protein = list(
#'     group1 = c("SHC1", "ERBB3", "ERBB2"),
#'     group2 = c("SRC", "EGFR")
#'   ),
#'   Drug = list(
#'     erlotinib = c("erlotinib")
#'   )
#' )
#'
#' # Plot Sankey diagram
#' GFusionMed::plot_sankey(edge_structure, edge_outcome, group_information)
#'
#' @export
plot_sankey <- function(edge_structure, edge_outcome, group_information, vector_layer_color = NULL) {

    if (!is.list(group_information)) {
        stop("group_information must be a list.")
    }

    vector_layer <- c(attr(edge_structure, "info")$structure_layer, attr(edge_outcome, "info")$outcome_layer)

    if (!is.null(vector_layer_color)) {
        
        if (!is.vector(vector_layer_color)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of structure layers.")
        }
        
        if (length(vector_layer_color) != length(vector_layer)) {
            stop("vector_layer_color must be a vector containing colors, and its length must match the total number of structure layers.")
        }
        
    }

    if (is.null(vector_layer_color)) {
        vector_layer_color <- assign_colors(vector_layer)
    }

    data_edge <- rbind(edge_structure, edge_outcome)

    data_edge <- subset(data_edge, Type == "dir")
    data_edge <- data_edge[, !(names(data_edge) %in% "Type")]

    data_edge <- transform(
        data_edge,
        Start_group = sapply(Start, assign_group, layer_info = group_information),
        End_group = sapply(End, assign_group, layer_info = group_information)
    )

    node_layer <-
        sapply(
            unique(c(data_edge$Start_group, data_edge$End_group)),
            function(x) strsplit(x, "_")[[1]][1]
        )

    node_colors <- setNames(vector_layer_color[node_layer], names(node_layer))

    node_vector <- names(node_layer)

    node_vector_index <- setNames(0:(length(node_vector)-1), node_vector)

    data_edge <- transform(
        data_edge,
        Start_group_index = sapply(Start_group, get_group_index, node_index = node_vector_index),
        End_group_index = sapply(End_group, get_group_index, node_index = node_vector_index),
        Color = paste0("rgba(128, 128, 128, ", Strength^2, ")")
    )

    fig <- 
        plotly::plot_ly(
            type = "sankey",
            orientation = "h",
            arrangement = "snap",
            
            node = list(
                label = node_vector,
                color = node_colors,
                pad = 15,
                thickness = 15,
                line = list(
                    color = "black",
                    width = 0.5
                )
            ),
            
            link = list(
                source = data_edge$Start_group_index,
                target = data_edge$End_group_index,
                value = rep(1, nrow(data_edge)),
                color = data_edge$Color
                )   
            )

    fig
}
