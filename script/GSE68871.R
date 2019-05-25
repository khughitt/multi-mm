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
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# exclude control sequences present in some datasets (GSE6477)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(sample_id=geo_accession, platform_id,
         response=`response to vtd therapy:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'BM-CD138+'

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = fData(eset)$`Gene Symbol`, .after = 1)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, sprintf('%s_expr.csv', accession)))
write_csv(sample_metadata, 
          file.path(clean_data_dir, sprintf('%s_sample_metadata.csv', accession)))

sessionInfo()

########################################################################################
#
# disabled..
#
# convert to gene symbols (for probes with multiple gene symbols, use the first one) 
# gene_symbols <- str_split(fData(eset)$`Gene Symbol`, '///', simplify=TRUE)
# gene_symbols <- str_trim(gene_symbols[, 1])
#
# convert to ensgene
#rownames(expr_dat) <- grch37$ensgene[match(gene_symbols, grch37$symbol)] 
#
#rownames(expr_dat) <- gene_symbols
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
