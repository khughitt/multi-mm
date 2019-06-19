#!/bin/env/Rscript
#
# Identification of multiple risk loci and regulatory mechanisms influencing susceptibility to multiple myeloma
# Went et al. (2018)
# n = 1038 samples
#
# - no survival, etc. related metadata fields.
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE31161'

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

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# report data processing used
print(as.character(pData(eset)$data_processing[1]))

# there are three samples in GSE31161 with all missing data..
num_missing <- apply(exprs(eset), 2, function(x) { sum(is.na(x)) })

#table(num_missing)
# num_missing
#     0 54675 
#  1035     3 

eset <- eset[, num_missing == 0]


# GSE31161 has not been adjusted for sample size
# range(colSums(exprs(eset), na.rm = TRUE))
# [1] 10358718 61443063

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

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

# unused
# treatment:ch1
# time of testing:ch1

# columns to include (GSE31161)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

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

