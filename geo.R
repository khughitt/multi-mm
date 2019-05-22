#!/bin/env/Rscript
library(annotables)
library(GEOquery)
library(tidyverse)

# geo accession
# https://www.ncbi.nlm.nih.gov/gds?Db=gds&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_gds&LinkReadableName=GEO+DataSets&ordinalpos=1&IdsFromResult=98717
# accession <- 'GSE6477'
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

# exclude control sequences present in some datasets (GSE6477)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# get relevant sample metadata
if (accession == 'GSE6477') {
  sample_metadata <- pData(eset) %>%
      select(sample_id=geo_accession, sample_title=title, ploidy=characteristics_ch1,
            ch13_status=characteristics_ch1.1) %>%
      mutate(mm_stage=ifelse(grepl('Normal', sample_title), 'Normal',
                      ifelse(grepl('New', sample_title), 'New',
                      ifelse(grepl('MGUS', sample_title), 'MGUS',
                      ifelse(grepl('Relapsed', sample_title), 'Relapsed',
                      ifelse(grepl('Smoldering', sample_title), 'Smoldering', 'Normal'))))))
} else if (accession == 'GSE83503') {
  # "treatment_protocol_ch1" is another alias for death;
  # "grow_protocol_ch1" is an alias for relapse;
  sample_metadata <- pData(eset) %>%
      select(sample_id=geo_accession, death=`death:ch1`, relapse=`relapse:ch1`)
}

# map to ensgene and collapse multi-gene mapped probes
if (accession == 'GSE6477') {
  # get expression data
  expr_dat <- exprs(eset)

  # convert to ensembl gene ids
  rownames(expr_dat) <- grch37$ensgene[match(fData(eset)$`Gene Symbol`, grch37$symbol)] 
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
} else if (accession == 'GSE83503') {
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
}


# store cleaned expression data and metadata
write.csv(expr_dat, file = file.path(clean_data_dir, sprintf('%s_expr.csv', accession)),
          quote = FALSE, row.names = FALSE)
write.csv(sample_metadata, file = file.path(clean_data_dir, 'sample_metadata.csv'),
          quote = FALSE, row.names = FALSE)

sessionInfo()
