#!/bin/env/Rscript
#
# Gene expression profiling of B lymphocytes and plasma cells from Waldenstrom's
# macroglobulinemia (2007)
#
# n = 52 patients 
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE6691'

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
# for GSE6691, data includes two separate esets with a small number of overlapping
# probes
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

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

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# columns to include (GSE6699)
# note: Healthy cells are from peripheral blood; all others are from bone marrow
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         disease_raw = characteristics_ch1)

sample_metadata$disease <- 'Healthy'
sample_metadata$disease[grepl('CLL', sample_metadata$disease_raw)] <- 'CLL'
sample_metadata$disease[grepl('MM', sample_metadata$disease_raw)] <- 'Multiple Myeloma'
sample_metadata$disease[grepl('WM', sample_metadata$disease_raw)] <- "Waldenstrom's macroglobulinemia"

sample_metadata <- sample_metadata %>%
  select(-disease_raw)

# add cell type and disease 
sample_metadata$cell_type = NA

# drop CLL samples (too different from MM)
mask <- sample_metadata$disease != 'CLL'

eset <- eset[, mask]
sample_metadata <- sample_metadata[mask, ]

# get gene symbols
gene_symbols <- fData(eset)$`Gene Symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

#table(expr_dat$gene_symbol == '')
# 
# FALSE  TRUE 
# 20599  1000 

# drop unmapped genes
expr_dat <- expr_dat %>%
  filter(gene_symbol != '')

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_expr.csv', accession)
sample_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))

sessionInfo()

