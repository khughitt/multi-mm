#!/bin/env/Rscript
library(GEOquery)
library(tidyverse)

options(stringsAsFactors = FALSE)

# GEO accession
accession <- 'GSE26760'

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

# GSE26760 has *not* been adjusted for sample size..
#range(colSums(exprs(eset)))
# [1] 6683420 9798074

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
#eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

# "Sample from patient MMRC0091" -> "MMRC0091"
patient_ids <- str_split(pData(eset)$source_name_ch1, ' ', simplify = TRUE)[, 4]

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

sample_metadata$patient_id <- patient_ids

# add cell type and disease (same for all samples)
sample_metadata$cell_type = 'CD138+'

# add additional metadata from
# http://portals.broadinstitute.org/mmgp/data/browseData?conversationPropagation=begin
mdat <- read_tsv(file.path(base_dir, 'metadata', 'mmrc.sample.information.txt')) %>%
	select(patient_id = Array, age = `Age at Diagnosis`, gender = Gender, 
         race = Race, mm_stage = Diagnosis)

sample_metadata <- sample_metadata %>%
	inner_join(mdat, by = 'patient_id') %>%
	filter(mm_stage != 'Unknown')

sample_metadata$mm_stage[grepl('Multiple Myeloma', sample_metadata$mm_stage)] <- 'MM'
sample_metadata$mm_stage[grepl('Leukemia', sample_metadata$mm_stage)] <- 'PCL'
sample_metadata$mm_stage[sample_metadata$mm_stage == 'Smoldering Myeloma'] <- 'SMM'

table(sample_metadata$mm_stage)
# 
# MGUS   MM  PCL  SMM 
#    2  224    3   10 
# 

# get gene symbols
gene_symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# drop samples with no available metadata
mask <- colnames(expr_dat) %in% c('probe_id', 'gene_symbol', sample_metadata$geo_accession)

#table(mask)
# mask
# FALSE  TRUE 
#    64   242 

expr_dat <- expr_dat[, mask]

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_expr.csv', accession)
sample_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, sample_outfile))

sessionInfo()

