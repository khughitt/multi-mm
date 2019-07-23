#!/bin/env/Rscript
#
# TODO July, 2019: 
#
# Find a probe -> gene mapping to use..
# BrainArray could potentially be used
#
#   https://www.ebi.ac.uk/arrayexpress/files/A-GEOD-4814/A-GEOD-4814_comments.txt
#
# However, there are no simple probe -> gene id mappings that I could find..
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE55145'

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
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

#
# data is missing...
#
# dim(eset)
# Features  Samples 
#        0       67 
#

# expression data must be downloaded separately
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55145/suppl/GSE55145_matrix_ifm67.txt.gz
expr_infile <- file.path(raw_data_dir, 'GSE55145_matrix_ifm67.txt.gz')

if (!file.exists(expr_infile)) {
  download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55145/suppl/GSE55145_matrix_ifm67.txt.gz',
                expr_infile)
}
expr_dat <- read.delim(gzfile(expr_infile), row.names = 1, skip = 4)

# report data processing used
#print(as.character(pData(eset)$data_processing[1]))

# [1] "HuEx-1_0-st-v2.r2.pgf"
# [1] "HuEx-1_0-st-v2.r2.pgf"

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

# columns to include (GSE55145)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         title,
         batch = `microarray_batch:ch1`,
         treatment_response = `post-induction (bortezomib) treatment response:ch1`) %>%
  mutate(sample_id = str_match(title, '^ifm67_\\d+')) %>%
  select(-title)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'CD138+'

# load expression data from file
#expr_dat <- read_tsv(expr_infile, skip = 4)
probe_ids <- expr_dat$u133pl2_id

expr_dat <- as.matrix(expr_dat[, grepl('ifm67_', colnames(expr_dat))])

#dim(fData(eset))
# [1]  0 12

rownames(expr_dat) <- probe_ids

# GSE55145 has *not* been adjusted for sample size..
range(colSums(expr_dat))
# [1] 249640.6 266477.4

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

# get gene symbols
gene_symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_expr.csv', accession)
sample_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))

sessionInfo()

