#!/bin/env/Rscript
library(GEOquery)
library(tidyverse)

options(stringsAsFactors = FALSE)

# GEO accession
accession <- 'GSE39754'

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

# GSE39754 has *not* been adjusted for sample size..
#range(colSums(exprs(eset)))
# [1] 2864927 4396094

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
#eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
#eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         diagnosis = `diagnosis:ch1`, treatment_response = `treatment_response:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$cell_type = 'CD138+'

sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$disease[grepl('Healthy', sample_metadata$diagnosis)] <- 'Healthy'

#all(sample_metadata$geo_accession == colnames(exprs(eset))) 
# [1] TRUE

# get gene symbols associated with each probe; gene symbols are stored at every
# [(N-1) + 2]th position (i.e. 2, 7, 12, 17..)
gene_parts <- str_split(fData(eset)$gene_assign, ' ///? ', simplify=TRUE)
gene_symbols <- gene_parts[, seq(2, ncol(gene_parts), by = 5)]

# collapse into the form "symbol 1 // symbol 2 // symbol 3 // ..."
gene_symbols <- apply(gene_symbols, 1, function(x) { str_trim(x)[x != ''] })
gene_symbols <- unlist(lapply(gene_symbols, paste, collapse = ' // '))

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# drop healthy doner samples
mask <- !sample_metadata$disease == 'Healthy'

sample_metadata <- sample_metadata[mask, ]
expr_dat <- expr_dat[, mask]

sample_metadata <- sample_metadata %>%
  select(-diagnosis)

# drop samples with no available metadata
#mask <- colnames(expr_dat) %in% c('probe_id', 'gene_symbol', sample_metadata$geo_accession)

#table(mask)
# mask
# TRUE 
#  178 

#expr_dat <- expr_dat[, mask]

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_expr.csv', accession)
sample_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))

sessionInfo()

