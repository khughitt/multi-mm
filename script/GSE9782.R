#!/bin/env/Rscript
#
# Gene expression profiling and correlation with outcome in clinical trials of the proteasome inhibitor bortezomib
# 
# Mulligan et al. (2007)
# n = 169 patients
#
# Description:
#
# The aims of this study were to assess the feasibility of prospective pharmacogenomics
# research in multicenter international clinical trials of bortezomib in multiple
# myeloma and to develop predictive classifiers of response and survival with
# bortezomib. Patients with relapsed myeloma enrolled in phase 2 and phase 3 clinical
# trials of bortezomib and consented to genomic analyses of pretreatment tumor samples.
# Bone marrow aspirates were subject to a negative-selection procedure to enrich for
# tumor cells, and these samples were used for gene expression profiling using DNA
# microarrays. Data quality and correlations with trial outcomes were assessed by
# multiple groups. Gene expression in this dataset was consistent with data published
# from a single-center study of newly diagnosed multiple myeloma. Response and survival
# classifiers were developed and shown to be significantly associated with outcome via
# testing on independent data. The survival classifier improved on the risk
# stratification provided by the International Staging System. Predictive models and
# biologic correlates of response show some specificity for bortezomib rather than
# dexamethasone. Informative gene expression data and genomic classifiers that predict
# clinical outcome can be derived from prospective clinical trials of new anticancer
# agents. Keywords: Gene expression profiling; correlation with outcome in clinical
# trials of the proteasome inhibitor bortezomib Overall design: Purified myeloma samples
# were collected prior to enrolment in clinical trials of bortezomib (PS-341). Samples
# were subject to replicate gene expression profiling using the Affymetrix 133A/B
# microarray. Data was normalized in MAS5.0 and the median of replicates is reported.
# Data was normalized to a Ttimmed mean of 15o and is NOT log transformed. Various
# patient parameters are reported as well as response, TTP and survival upon treatment
# with bortezomib or dexamethasone.
# 
# Source: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA103715
#
library(annotables)
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE9782'

# directory to store raw and processed data
base_dir <- file.path('/data/public/human/geo', accession)

raw_data_dir <- file.path(base_dir, 'raw')
clean_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, clean_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
# for GSE9782, data includes two separate esets with a small number of overlapping
# probes
esets <- getGEO(accession, destdir = raw_data_dir)

for (i in 1:length(esets)) {
  eset <- esets[[i]]

  # list the columns that are neither all different or all unique
  covariates <- c()

  for (cname in colnames(pData(eset))) {
    num_vals <- length(table(pData(eset)[, cname]))

    if (num_vals != 1 && num_vals != ncol(eset)) {
      covariates <- c(covariates, cname)

      message(cname)
      print(table(pData(eset)[, cname]))
    }
  }

  # exclude control sequences present in some datasets
  eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

  # columns to include (GSE9782)
  # characteristics_ch1.1 (treatment)
  # characteristics_ch1.6 (myeloma measure)
  # characteristics_ch1.7 (pgx_response)
  # characteristics_ch1.8 (pgx_responder)

  # unused columns
  #
  # characteristics_ch1 (study code)
  # characteristics_ch1.2 (gender)
  # characteristics_ch1.3 (ethnicity)
  #

  # 
  # other potentially interesting information that is not consistently encoded:
  #
  # characteristics_ch1.7     <fct> PGx_Response = NC                                          
  # characteristics_ch1.8     <fct> PGx_Responder = NR                                         
  # characteristics_ch1.9     <fct> PGx_Days_To_Progression = 87                               
  # characteristics_ch1.10    <fct> PGx_CensorReason = No more data                            
  # characteristics_ch1.11    <fct> "Did_Patient_Die(0=No,1=Yes) = 1"                          
  # characteristics_ch1.12    <fct> Days_Survived_From_Randomization = 106                     
  # characteristics_ch1.13    <fct> albumin = 41                              
  #
  # > pData(eset)[, 'characteristics_ch1.9']                                                     
  # [1] PGx_Days_To_Progression = 87           PGx_Progression(0=No,1=Yes) = 1                 
  # [3] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1                 
  # [5] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1                 
  # [7] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1                 
  # [9] PGx_CensorReason = No more data        PGx_Progression(0=No,1=Yes) = 1  

  # get relevant sample metadata
  # characteristics_ch1.1 (treatment)
  # characteristics_ch1.6 (myeloma measure)
  # characteristics_ch1.7 (pgx_response)
  # characteristics_ch1.8 (pgx_responder)
  sample_metadata <- pData(eset) %>%
    select(sample_id = geo_accession, platform_id, treatment = characteristics_ch1.1,
           myeloma_measure = characteristics_ch1.6, pgx_response = characteristics_ch1.7,
           pgx_responder = characteristics_ch1.8)

  # add cell type and disease (same for all samples)
  sample_metadata$disease = 'Multiple Myeloma'
  sample_metadata$cell_type = NA # likely CD138+, but not 100% clear

  # get gene symbols and replace '///' with '//' to be consistent with other datasets
  gene_symbols <- str_replace(fData(eset)$`Gene Symbol`, '///', '//') 

  # get expression data and add gene symbol column
  expr_dat <- exprs(eset) %>%
    as.data.frame %>%
    rownames_to_column('probe_id') %>%
    add_column(gene_symbol = gene_symbols, .after = 1)

  # determine filenames to use for outputs and save to disk
  if (i == 1) {
    expr_outfile <- sprintf('%s_expr.csv', accession)
    sample_outfile <- sprintf('%s_sample_metadata.csv', accession)
  } else {
    expr_outfile <- sprintf('%s_expr_%d.csv', accession, i)
    sample_outfile <- sprintf('%s_sample_metadata_%d.csv', accession, i)
  }

  # store cleaned expression data and metadata
  write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
  write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))
}

sessionInfo()

########################################################################################
#
# disabled..
#
# convert to ensembl gene ids (use first symbol listed..)
#gene_symbols <- str_split(fData(eset)$`Gene Symbol`, '///', simplify=TRUE)
#gene_symbols <- str_trim(gene_symbols[, 1])
#
#rownames(expr_dat) <- grch37$ensgene[match(gene_symbols, grch37$symbol)] 
#
#rownames(expr_dat) <- gene_symbols
#
#mask <- !is.na(rownames(expr_dat))
#message(sprintf("Dropping %d / %d probes which could not be mapped to gene symbols...",
#                sum(!mask), nrow(expr_dat)))
#expr_dat <- expr_dat[mask, ]
#
# sum multi-mapped gene ids
#num_before <- nrow(expr_dat)
#expr_dat <- aggregate(expr_dat, list(rownames(expr_dat)), sum)
#rownames(expr_dat) <- expr_dat[, 1]
#expr_dat <- expr_dat[, -1]
#message(sprintf("%d / %d genes remain after averaging.", nrow(expr_dat), num_before))
#
########################################################################################
