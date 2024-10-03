
# -*- coding: utf-8 -*-

rm(list = ls())

# install.packages("BiocManager")

# BiocManager::install("depmap")

suppressPackageStartupMessages({
    library(depmap)
})

path_data <- "/root/Project/GFusionMed/application/data"

# CNA (depmap 22Q2)
cna <- depmap::depmap_copyNumber()

# mRNA (depmap 22Q2)
mrna <- depmap::depmap_TPM()

# protein (depmap 19Q3)
protein <- depmap::depmap_RPPA()

# Meta (depmap 22Q2)
meta <- depmap::depmap_metadata()

# Drug (depmap 21Q2)
drug <- depmap::depmap_drug_sensitivity()

env_save <- new.env(parent = emptyenv())
env_save$cna <- cna
env_save$mrna <- mrna
env_save$protein <- protein
env_save$meta <- meta
env_save$drug <- drug

do.call(
    "save",
    c(
        ls(envir = env_save),
        list(envir = env_save, file = paste0(path_data, "/data_cna_mrna_protein_meta_drug.RData"))
    )
)
