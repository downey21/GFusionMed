
# -*- coding: utf-8 -*-

rm(list = ls())

suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(tibble)
})

path_data <- "/root/Project/GFusionMed/application/data"

loaded_data <- load(paste0(path_data, "/data_cna_mrna_protein_meta_drug.RData"))

ab_map <- tibble::as_tibble(read.csv(paste0(path_data, "/CCLE_RPPA_Ab_info_20181226.csv")))

drug_lung_cancer <-
    c(
        "erlotinib",
        "gefitinib",
        "afatinib",
        "dacomitinib",
        "osimertinib",
        "cisplatin"
    )

# common depmap_id
lung_cancer_depmap_id <-
    meta %>%
    dplyr::filter(primary_disease == "Lung Cancer") %>%
    dplyr::pull(depmap_id)

lung_cancer_depmap_id_cna <-
    cna %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id) %>%
    dplyr::pull(depmap_id)

lung_cancer_depmap_id_mrna <-
    mrna %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id) %>%
    dplyr::pull(depmap_id)

lung_cancer_depmap_id_protein <-
    protein %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id) %>%
    dplyr::pull(depmap_id)

lung_cancer_depmap_id_drug <-
    drug %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id) %>%
    dplyr::pull(depmap_id)

lung_cancer_depmap_id_intersect <-
    base::intersect(
        base::intersect(
            base::intersect(lung_cancer_depmap_id_cna, lung_cancer_depmap_id_mrna),
            lung_cancer_depmap_id_protein
        ),
        lung_cancer_depmap_id_drug
    )

# make dataset for lung cancer data
cna_lung_cancer <-
    cna %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id_intersect) %>%
    dplyr::select(c("depmap_id", "gene_name", "log_copy_number")) %>% 
    tidyr::pivot_wider(names_from = gene_name, values_from = log_copy_number) %>%
    dplyr::arrange(depmap_id)

mrna_lung_cancer <-
    mrna %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id_intersect) %>%
    dplyr::select(c("depmap_id", "gene_name", "rna_expression")) %>% 
    tidyr::pivot_wider(names_from = gene_name, values_from = rna_expression) %>%
    dplyr::arrange(depmap_id)

protein_lung_cancer <-
    protein %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id_intersect) %>%
    dplyr::select(c("depmap_id", "antibody", "expression")) %>% 
    tidyr::pivot_wider(names_from = antibody, values_from = expression) %>%
    dplyr::arrange(depmap_id)

colnames(protein_lung_cancer) <- c("depmap_id", unname(sapply(ab_map$Target_Genes, function(x) stringr::str_split(x, " ")[[1]][1])))
protein_lung_cancer <- protein_lung_cancer[, !duplicated(names(protein_lung_cancer))]

drug_lung_cancer <-
    drug %>%
    dplyr::filter(depmap_id %in% lung_cancer_depmap_id_intersect) %>%
    dplyr::select(c("depmap_id", "name", "dependency")) %>%
    dplyr::filter(name %in% drug_lung_cancer) %>% 
    tidyr::pivot_wider(names_from = name, values_from = dependency) %>%
    dplyr::arrange(depmap_id)

head(cna_lung_cancer)
head(mrna_lung_cancer)
head(protein_lung_cancer)
head(drug_lung_cancer)

env_save <- new.env(parent = emptyenv())
env_save$cna_lung_cancer <- cna_lung_cancer
env_save$mrna_lung_cancer <- mrna_lung_cancer
env_save$protein_lung_cancer <- protein_lung_cancer
env_save$drug_lung_cancer <- drug_lung_cancer

do.call(
    "save",
    c(
        ls(envir = env_save),
        list(envir = env_save, file = paste0(path_data, "/data_cna_mrna_protein_drug_lung_cancer.RData"))
    )
)
