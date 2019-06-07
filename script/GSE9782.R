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
esets <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)

# report data processing used
print(as.character(pData(esets[[1]])$data_processing[1]))

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

  # GSE9782 has *not* been adjusted for sample size..
  #range(colSums(exprs(eset)))
  # [1] 5772391 9324047

  # perform size-factor normalization
  exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

  # exclude control sequences present in some datasets
  eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

  # 
  # Note: some of the metadata fields appear to be incorrectly encoded, e.g.:
  #
  # > pData(eset)[, 'characteristics_ch1.9']                                                     
  # [1] PGx_Days_To_Progression = 87           PGx_Progression(0=No,1=Yes) = 1                 
  # [3] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1                 
  # [5] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1                 
  # [7] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1                 
  # [9] PGx_CensorReason = No more data        PGx_Progression(0=No,1=Yes) = 1  

  #
  # to get around this, we will manually detect and parse those fields..
  #

  sample_metadata <- pData(eset) %>%
    mutate(
      study_code      = as.numeric(str_match(characteristics_ch1, '\\d+$')),
      treatment       = str_match(characteristics_ch1.1, '\\w+$'),
      gender          = str_match(characteristics_ch1.2, '\\w+$'),
      ethnicity       = str_match(characteristics_ch1.3, '\\w+$'),
      age             = as.numeric(str_match(characteristics_ch1.4, '\\d+$')),
      pgx_response    = str_match(characteristics_ch1.7, '\\w+$'), 
      pgx_responder   = str_match(characteristics_ch1.8, '\\w+$')) %>%
    select(sample_id = geo_accession, platform_id, study_code, treatment, gender,
           ethnicity, age, pgx_response, pgx_responder)

  # convert metadata from factor to character columns for parsing
  str_mdat <- pData(eset)

  for (cname in colnames(str_mdat)) {
    str_mdat[, cname] <- as.character(str_mdat[, cname])
  }

  pfs_censor <- c()
  pfs_days <- c()
  pfs_censor_reason <- c()
  os_censor <- c()
  os_days <- c()

  # iterate over samples and find relevant information
  for (sample_num in 1:nrow(pData(eset))) {
    # PGx Progression
    ind <- which(grepl('PGx_Prog', str_mdat[sample_num, ]))

    if (length(ind) == 0) {
      pfs_censor <- c(pfs_censor, 0)
    } else {
      pfs_censor <- c(pfs_censor, ifelse(endsWith(str_mdat[sample_num, ind], '1'), 1, 0))
    }

    # PGx days
    ind <- which(grepl('PGx_Days', str_mdat[sample_num, ]))

    if (length(ind) == 0) {
      pfs_days <- c(pfs_days, NA)
    } else {
      pfs_days <- c(pfs_days, as.numeric(str_match(str_mdat[sample_num, ind], '\\d+$')))
    }

    # PGx Censor Reason
    ind <- which(grepl('PGx_CensorReason', str_mdat[sample_num, ]))

    if (length(ind) == 0) {
        pfs_censor_reason <- c(pfs_censor_reason, NA)
    } else {
      pfs_censor_reason <- c(pfs_censor_reason, trimws(str_match(str_mdat[sample_num, ind], '[ \\w]+$')))
    }

    # Deceased
    ind <- which(grepl('Did_Patient_Die', str_mdat[sample_num, ]))

    if (length(ind) == 0) {
      os_censor <- c(os_censor, 0)
    } else {
      os_censor <- c(os_censor, 1)
    }

    # Days Survived
    ind <- which(grepl('Days_Survived', str_mdat[sample_num, ]))
    os_days <- c(os_days, as.numeric(str_match(str_mdat[sample_num, ind], '\\d+$')))
  }

  # add to sample metadata
  sample_metadata$pfs_days <- pfs_days
  sample_metadata$pfs_censor <- pfs_censor
  sample_metadata$pfs_censor_reason <- pfs_censor_reason
  sample_metadata$os_days <- os_days
  sample_metadata$os_censor <- os_censor

  # add cell type and disease (same for all samples)
  sample_metadata$disease = 'Multiple Myeloma'
  sample_metadata$cell_type = NA # likely CD138+, but not 100% clear

  # get gene symbols
  gene_symbols <- fData(eset)$`Gene symbol`

  # get expression data and add gene symbol column
  expr_dat <- exprs(eset) %>%
    as.data.frame %>%
    rownames_to_column('probe_id') %>%
    add_column(gene_symbol = gene_symbols, .after = 1)

  # determine filenames to use for outputs and save to disk
  expr_outfile <- sprintf('%s_%d_expr.csv', accession, i)
  sample_outfile <- sprintf('%s_%d_sample_metadata.csv', accession, i)

  # store cleaned expression data and metadata
  write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
  write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))
}

sessionInfo()

