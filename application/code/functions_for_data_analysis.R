
load_data_for_structure <- function(path_data, pp) {
        loaded_data <- load(paste0(path_data,"/data_cna_mrna_protein_drug_lung_cancer.RData"))

        CNA <- cna_lung_cancer
        mRNA <- mrna_lung_cancer
        RPPA <- protein_lung_cancer

        load(paste0(path_data,"/StringV10_ppi_pathway1.RData"))
        padat <- ppi$pathwaydat
        padat[,2] <- gsub("-","",padat[,2])
        pname <- unique(padat[,1])[pp]
        genes <- unlist(strsplit(padat[padat[,1] == pname,3],split=", "))

        RPPA <- RPPA[,colnames(RPPA) %in% genes]
        mRNA <- mRNA[,colnames(mRNA) %in% genes]
        CNA <- CNA[,colnames(CNA) %in% genes]

        data_ready <-
            list(
                CNA = CNA,
                mRNA = mRNA,
                Protein = RPPA
            )

        return(data_ready)
    }

load_data_for_outcome <- function(path_data, pp, drug_name) {
        loaded_data <- load(paste0(path_data,"/data_cna_mrna_protein_drug_lung_cancer.RData"))

        CNA <- cna_lung_cancer
        mRNA <- mrna_lung_cancer
        RPPA <- protein_lung_cancer
        Drug <- drug_lung_cancer

        Drug <- Drug[,c("depmap_id", drug_name)]
        Drug <- Drug[rowSums(is.na(Drug))==0,]

        rnames <- Drug$depmap_id
        CNA <- CNA[match(rnames,CNA$depmap_id),]
        mRNA <- mRNA[match(rnames,mRNA$depmap_id),]
        RPPA <- RPPA[match(rnames,RPPA$depmap_id),]

        load(paste0(path_data,"/StringV10_ppi_pathway1.RData"))
        padat <- ppi$pathwaydat
        padat[,2] <- gsub("-","",padat[,2])
        pname <- unique(padat[,1])[pp]
        genes <- unlist(strsplit(padat[padat[,1] == pname,3],split=", "))

        RPPA <- RPPA[,colnames(RPPA) %in% genes]
        mRNA <- mRNA[,colnames(mRNA) %in% genes]
        CNA <- CNA[,colnames(CNA) %in% genes]
        Drug <- Drug[-1]

        data_ready <-
            list(
                CNA = CNA,
                mRNA = mRNA,
                Protein = RPPA,
                Drug = Drug
            )

        return(data_ready)
    }

mediation_summary_all_for_each_outcome <- function(path_data, path_result, drug_name) {
    
    load(paste0(path_data,"/StringV10_ppi_pathway1.RData"))
    padat = ppi$pathwaydat
    padat[,2] = gsub("-","",padat[,2])
    ppname = unique(padat[,1])

    result <- NULL
    for (pp in 1:length(ppname)) {

        load(paste(path_result, "/fit_structure_pp", pp, ".RData", sep = ""))
        result_structure <- result_fit
        
        load(paste(path_result,"/fit_outcome_", drug_name, "_pp", pp, ".RData", sep = ""))
        result_outcome <- result_fit

        result_temp <- GFusionMed::perform_mediation_analysis(result_structure, result_outcome)

        data_frame_temp <- data.frame(matrix(ppname[pp], nrow = nrow(result_temp), ncol = 1))
        colnames(data_frame_temp) <- "Pathway"
        
        result_temp <- cbind(data_frame_temp, result_temp)

        result <- rbind(result, result_temp)
    }

    return(result)
}

mediation_summary_all <- function(path_data, path_result, vector_drug) {
    result <- NULL

    for (drug_name in vector_drug) {
        result_temp <- mediation_summary_all_for_each_outcome(path_data, path_result, drug_name)

        result <- rbind(result, result_temp)
    }

    return(result)
}

edge_summary_all_for_each_outcome <- function(path_data, path_result, drug_name) {
    
    load(paste0(path_data,"/StringV10_ppi_pathway1.RData"))
    padat = ppi$pathwaydat
    padat[,2] = gsub("-","",padat[,2])
    ppname = unique(padat[,1])

    result <- NULL
    for (pp in 1:length(ppname)) {

        load(paste(path_result, "/fit_structure_pp", pp, ".RData", sep = ""))
        result_structure <- result_fit
        
        load(paste(path_result,"/fit_outcome_", drug_name, "_pp", pp, ".RData", sep = ""))
        result_outcome <- result_fit

        result_temp_structure <- GFusionMed::edge_summary(result_structure)
        result_temp_outcome <- GFusionMed::edge_summary(result_outcome)

        result_temp <- rbind(result_temp_structure, result_temp_outcome)

        data_frame_temp <- data.frame(matrix(ppname[pp], nrow = nrow(result_temp), ncol = 1))
        colnames(data_frame_temp) <- "Pathway"
        
        result_temp <- cbind(data_frame_temp, result_temp)

        result <- rbind(result, result_temp)
    }

    return(result)
}

edge_summary_all <- function(path_data, path_result, vector_drug) {
    result <- NULL

    for (drug_name in vector_drug) {
        result_temp <- edge_summary_all_for_each_outcome(path_data, path_result, drug_name)

        data_frame_temp <- data.frame(matrix(drug_name, nrow = nrow(result_temp), ncol = 1))
        colnames(data_frame_temp) <- "Drug"

        result_temp <- cbind(data_frame_temp, result_temp)

        result <- rbind(result, result_temp)
    }

    return(result)
}
