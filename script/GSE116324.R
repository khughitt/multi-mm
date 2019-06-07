#!/bin/env/Rscript
#
# RNA-Seq of newly diagnosed patients in the PADIMAC study leads to a bortezomib/lenalidomide decision signature
# Chapman et al. (2018)
# n = 44
#
library(annotables)
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE116324'

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
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

#dim(eset)
# Features  Samples 
#        0       44 

# getGEO only retrieves sample metadata, but no expression data or gene metadata;
# counts are available as a separate supplementary file, however..
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116324/suppl/GSE116324_padimacRnaSeq.csv.gz
supp_file <- file.path(raw_data_dir, accession, 'GSE116324_padimacRnaSeq.csv.gz')

if (!file.exists(supp_file)) {
  getGEOSuppFiles(accession, baseDir = raw_data_dir)[[1]]
}
expr <- read.csv(gzfile(supp_file), row.names = 1)

#expr[1:3, 1:3]
#                   PAD.004 PAD.017 PAD.018
# ENSG00000223972.4       0       0       0
# ENSG00000227232.4      45    1242     223
# ENSG00000243485.2       0       0       0

rownames(expr) <- sub('\\.\\d+', '', rownames(expr)) 

# report data processing used
print(as.character(pData(eset)$data_processing[1]))

# GSE116324 has not been adjusted for sample size
#range(colSums(expr))
# [1]  2912575 14114474

# perform size-factor normalization
expr <- sweep(expr, 2, colSums(expr), '/') * 1E6

# exclude zero variance genes
row_vars <- apply(expr, 1, var)

#sum(row_vars == 0)
# [1] 18030

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


# unused
# treatment:ch1
# time of testing:ch1

# columns to include (GSE116324)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
  bortezomib_response = `bortezomib response:ch1`,
  iss_stage = `iss stage:ch1`,
  age = `age:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'CD138+'

# get gene symbols
gene_symbols <- grch37$symbol[match(rownames(expr), grch37$ensgene)]

# get expression data and add gene symbol column
expr_dat <- expr %>%
  rownames_to_column('ensgene') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_expr.csv', accession)
sample_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))

sessionInfo()

