#!/bin/env/Rscript
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE47552'

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

# GSE47552 has *not* been adjusted for sample size..
#range(colSums(exprs(eset)))
# [1] 41046258 65576570

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
#eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
#eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

# columns to include (GSE47552)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         mm_stage_raw = `cell type:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'CD138+'

sample_metadata$disease[grepl('NPC', sample_metadata$mm_stage_raw)] <- 'Healthy'

sample_metadata$mm_stage <- rep('MM', length(sample_metadata$mm_stage_raw))

sample_metadata$mm_stage[grepl('Normal', sample_metadata$mm_stage_raw)] <- 'Healthy'
sample_metadata$mm_stage[grepl('MGUS', sample_metadata$mm_stage_raw)] <- 'MGUS'
sample_metadata$mm_stage[grepl('SMM', sample_metadata$mm_stage_raw)] <- 'SMM'

sample_metadata <- sample_metadata %>%
  select(-mm_stage_raw)

# get gene symbols
gene_symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# table(expr_dat$gene_symbol == '')
# 
# FALSE  TRUE 
# 22195 11102 

# drop unmapped genes
expr_dat <- expr_dat %>%
  filter(gene_symbol != '')

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.csv', accession)
mdat_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, mdat_outfile))

sessionInfo()

