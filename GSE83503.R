#!/bin/env/Rscript
library(annotables)
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE83503'

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

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(sample_id=geo_accession, platform_id, death=`death:ch1`, relapse=`relapse:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'BM-CD138+'

# map to ensgene and collapse multi-gene mapped probes

# exclude control spots (4130 / 22011 samples)
mask <- fData(eset)$SPOT_ID != 'control'
eset <- eset[mask, ]

# convert to gene symboles
# "NM_001130045 // TTLL10 // tubulin tyrosine ligase-like family, member 10 // ..."
GENE_SYMBOL_IDX <- 2
gene_symbols <- str_split(fData(eset)$gene_assignment, '//', simplify=TRUE)
gene_symbols <- str_trim(gene_symbols[, GENE_SYMBOL_IDX])

# get expression matrix
expr_dat <- exprs(eset)

# map from gene symbols -> ensgene
rownames(expr_dat) <- grch37$ensgene[match(gene_symbols, grch37$symbol)]

# drop any unmappable ids
mask <- !is.na(rownames(expr_dat))
message(sprintf("Dropping %d / %d probes which could not be mapped to Ensembl genes...",
                sum(!mask), nrow(expr_dat)))
expr_dat <- expr_dat[mask, ]

# store cleaned expression data and metadata
write.csv(expr_dat, file = file.path(clean_data_dir, sprintf('%s_expr.csv', accession)),
          quote = FALSE, row.names = FALSE)
write.csv(sample_metadata, file = file.path(clean_data_dir, 'sample_metadata.csv'),
          quote = FALSE, row.names = FALSE)

sessionInfo()
