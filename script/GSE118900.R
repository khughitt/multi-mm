#!/bin/env/Rscript
#
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE118900'

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
#
# eset for this dataset only include metadata; expression data is empty and will be
# loaded separately..
#
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

expr <- read.delim(gzfile('/data/public/human/geo/GSE118900/raw/GSE118900_MM.scrna-seq.tpm.pass.txt.gz'), row.names=1)

# replace rownames with gene symbol

#head(rownames(expr))
# [1] "AADACL3|chr1|12776118" "AADACL4|chr1|12704566" "ABCA4|chr1|94458394"   "ABCB10|chr1|229652329"
# [5] "ABCD3|chr1|94883933"   "ABL2|chr1|179068462"  

symbols <- str_split(rownames(expr), '\\|', simplify = TRUE)[, 1]

# remove small number of duplicated gene symbol entries
mask <- !duplicated(symbols)

expr <- expr[mask, ]
rownames(expr) <- symbols[mask]

# report data processing used
print(as.character(pData(eset)$data_processing[1]))

# GSE118900 has been adjusted for sample size
#range(colSums(expr))
# [1]  993252.9 1000000.0

# perform size-factor normalization
expr <- sweep(expr, 2, colSums(expr), '/') * 1E6

# exclude zero variance genes
row_vars <- apply(expr, 1, var)

#sum(row_vars == 0)
# [1] 3760

expr <- expr[row_vars != 0, ]

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

# columns to include (GSE118900)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         patient = `patient:ch1`,
         mm_stage = `tumor stage:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'CD138+'

expr <- expr %>%
  rownames_to_column('gene_symbol')

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.csv', accession)
mdat_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, mdat_outfile))

sessionInfo()

