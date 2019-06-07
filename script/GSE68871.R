#!/bin/env/Rscript
library(annotables)
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE68871'

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
# result is a list with a single entry containing an ExpressionSet instance
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# report data processing used
print(pData(eset)$data_processing[1])

# GSE68871 has already been size-factor scaled
#range(colSums(exprs(eset)))
# [1] 288678.1 292655.8

# in order to normalize downstream comparisons across datasets, however, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets (GSE6477)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# GSE68871 has no zero-variance probes..
#any(apply(exprs(eset), 1, var) == 0)
# [1] FALSE

# exclude any probes with zero variance (uninformative)
#eset <- eset[apply(exprs(eset), 1, var) > 0, ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         treatment_response=`response to vtd therapy:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'BM-CD138+'

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = fData(eset)$`Gene symbol`, .after = 1)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, sprintf('%s_expr.csv', accession)))
write_csv(sample_metadata, 
          file.path(clean_data_dir, sprintf('%s_sample_metadata.csv', accession)))

sessionInfo()
