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

# map to ensgene and collapse multi-gene mapped probes
expr_dat <- exprs(eset)

# convert to ensembl gene ids (use first symbol listed..)
gene_symbols <- str_split(fData(eset)$`Gene Symbol`, '///', simplify=TRUE)
gene_symbols <- str_trim(gene_symbols[, 1])

rownames(expr_dat) <- grch37$ensgene[match(gene_symbols, grch37$symbol)] 
mask <- !is.na(rownames(expr_dat))
message(sprintf("Dropping %d / %d probes which could not be mapped to Ensembl genes...",
                sum(!mask), nrow(expr_dat)))
expr_dat <- expr_dat[mask, ]

# sum multi-mapped gene ids
num_before <- nrow(expr_dat)
expr_dat <- aggregate(expr_dat, list(rownames(expr_dat)), sum)
rownames(expr_dat) <- expr_dat[, 1]
expr_dat <- expr_dat[, -1]
message(sprintf("%d / %d genes remain after averaging.", nrow(expr_dat), num_before))

# store cleaned expression data and metadata
write.csv(expr_dat, file = file.path(clean_data_dir, sprintf('%s_expr.csv', accession)),
          quote = FALSE, row.names = FALSE)
write.csv(sample_metadata, file = file.path(clean_data_dir, 'sample_metadata.csv'),
          quote = FALSE, row.names = FALSE)

sessionInfo()
