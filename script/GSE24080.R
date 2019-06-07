#!/bin/env/Rscript
#
# MAQC-II Project: Multiple myeloma (MM) data set (GSE24080)
# Popovici et al. (2010)
# n = 559 patients
#
# Note: this dataset appears to use the same samples from GSE2658, but processed in
# a different manner, and with different metadata.
#
# Interestingly, the MAQC-II version of the dataset includes ~2x samples after filtering
# and also includes survival-related metadata.
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE24080'

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
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# report data processing used
print(as.character(pData(eset)$data_processing[1]))

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

# GSE24080 has *not* been adjusted for sample size..
#range(colSums(exprs(eset)))
# [1] 367378.3 466299.2

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

# exclude samples with no expression data
num_missing <- apply(exprs(eset), 2, function(x) { sum(is.na(x)) })

#table(num_missing)
# num_missing
#   0 
# 559 

# unused

# columns to include (GSE24080)
efs <- as.numeric(endsWith(pData(eset)[, 'efs milestone outcome (24 months):ch1'], '0')) 
os  <- as.numeric(endsWith(pData(eset)[, 'os milestone outcome (24 months):ch1'], '0')) 

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         gender=`Sex:ch1`,
         age=`age:ch1`,
         maqc_status=`maqc_distribution_status:ch1`) %>%
  add_column(efs_event = efs, os_event = os, .after = 'platform_id')

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'CD138+'

# get gene symbols
gene_symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_expr.csv', accession)
sample_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))

sessionInfo()

