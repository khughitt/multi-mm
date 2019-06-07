#!/bin/env/Rscript
#
# Gene Expression Profiles of Multiple Myeloma (GSE2658) /
# MAQC-II Project: Multiple myeloma (MM) data set (GSE24080)
# 
# Hanamura et al. (2006)
# Popovici et al. (2010)
#
# n = 559 patients
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE2658'
#accession <- 'GSE24080'

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

# GSE2658 has *not* been adjusted for sample size..
#range(colSums(exprs(eset)))
# [1] 5772391 9324047

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
#     0 54613 
#   378   181 

eset <- eset[, num_missing == 0]

# columns to include (GSE2658)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         deceased = characteristics_ch1,
         patient_subgroup = characteristics_ch1.8, characteristics_ch1.10) %>%
         mutate(contaminated = characteristics_ch1.10 == '[Subgrp7=MY]') %>%
         mutate(deceased = startsWith(as.character(deceased), '[SURIND=1')) %>%
         select(-characteristics_ch1.10) %>%
         mutate(patient_subgroup = sub("Subgrp7=", "", str_extract(patient_subgroup, "Subgrp7=[[:alnum:]]+")))

# exclude contaminated samples
sample_metadata <- sample_metadata %>%
  filter(!contaminated) %>%
  select(-contaminated)

mask <- colnames(eset) %in% sample_metadata$geo_accession

#table(mask)
# mask
# FALSE  TRUE 
#    95   283 

eset <- eset[, mask]

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

