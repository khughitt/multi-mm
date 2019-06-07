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

# report data processing used
print(as.character(pData(eset)$data_processing[1]))

# exclude control spots (4130 / 22011 samples)
mask <- fData(eset)$SPOT_ID != 'control'
eset <- eset[mask, ]

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var) > 0, ]

# GSE83503 has already been size-factor scaled, so no need to adjust
#range(colSums(exprs(eset)))
# [1]  97004.52 106777.29

# in order to normalize downstream comparisons across datasets, however, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# GSE83503 has no zero-variance probes..
# any(apply(exprs(eset), 1, var) == 0)
# [1] FALSE

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, death=`death:ch1`, relapse=`relapse:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'BM-CD138+'

# get gene symbols associated with each probe; gene symbols are stored at ever
# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fData(eset)$gene_assignment, '//', simplify=TRUE)
gene_symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
gene_symbols <- apply(gene_symbols, 1, function(x) { str_trim(x)[x != ''] })
gene_symbols <- unlist(lapply(gene_symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, sprintf('%s_expr.csv', accession)))
write_csv(sample_metadata, 
          file.path(clean_data_dir, sprintf('%s_sample_metadata.csv', accession)))

sessionInfo()

########################################################################################
#
# disabled..
#
# convert to gene symboles
# "NM_001130045 // TTLL10 // tubulin tyrosine ligase-like family, member 10 // ..."
#GENE_SYMBOL_IDX <- 2
#gene_symbols <- str_split(fData(eset)$gene_assignment, '//', simplify=TRUE)
#gene_symbols <- str_trim(gene_symbols[, GENE_SYMBOL_IDX])
#
## map from gene symbols -> ensgene
##rownames(expr_dat) <- grch37$ensgene[match(gene_symbols, grch37$symbol)]
#
#rownames(expr_dat) <- gene_symbols
#
## drop any unmappable ids
#mask <- !is.na(rownames(expr_dat))
#message(sprintf("Dropping %d / %d probes which could not be mapped to gene symbols...",
#                sum(!mask), nrow(expr_dat)))
#expr_dat <- expr_dat[mask, ]
#
########################################################################################
