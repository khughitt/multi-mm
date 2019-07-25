#!/bin/env/Rscript
#
# Gene expression profiling of B lymphocytes and plasma cells from Waldenstrom's
# macroglobulinemia (2007)
#
# n = 52 patients 
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE6691'

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

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

#head(colnames(eset))
# [1] "5.82592"  "6.799907" "5.388776" "6.630125" "7.093243" "6.426071"

# NOTE: GEOquery currently incorrectly parses the data for this experiment, skipping
# >500 lines and usign the wrong colnames; reported the issue upstream. For now, just
# fixing colnames..
colnames(eset) <- pData(eset)$geo_accession

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

# columns to include (GSE6699)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

# add cell type and disease 
sample_metadata$cell_type = NA

#table(pData(eset)$characteristics_ch1)
# 
#                                                                                         BL from CLL patient, isolated using CD19-PE and CD5-APC 
#                                                                                                                                               2 
#                                                                                        BL from CLL patient, isolated using CD19-PE and CD5-APC. 
#                                                                                                                                               9 
#                                                                                             BL from healthy donors, isolated using CD19-PE-Cy7. 
#                                                                                                                                               8 
# BL from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE, CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7. 
#                                                                                                                                               9 
#    from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE, CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7. 
#                                                                                                                                               1 
#                                                                                                PC from healthy donors, isolated using CD38-APC. 
#                                                                                                                                               5 
#                                                                                                     PC from MM patient, isolated using CD38-APC 
#                                                                                                                                               2 
#                                                                                                    PC from MM patient, isolated using CD38-APC. 
#                                                                                                                                               9 
#                                                                                                     PC from MM patient,isolated using CD38-APC. 
#                                                                                                                                               1 
# PC from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE, CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7. 
#                                                                                                                                               9 
#  PC from WM patient, isolated using Kappa or Lambda-fluorescein isothiocyanate, CD10-PE,CD38-PerCPCy5.5,CD19-PE-Cy7, CD34-APC and CD45-APC-Cy7. 
#                                                                                                                                               1 
desc <- pData(eset)$characteristics_ch1

disease <- rep('Multiple Myeloma', length(desc))
disease[grepl('healthy', desc)] <- "Healthy"
disease[grepl('WM', desc)] <- "Waldenström's Macroglobulinemia"
disease[grepl('CLL', desc)] <- "Chronic Lymphocytic Leukemia"

#table(disease)
# disease
#    Chronic Lymphocytic Leukemia                         Healthy                Multiple Myeloma 
#                              11                              13                              12 
# Waldenström's Macroglobulinemia 
#                              20 
sample_metadata$disease <- disease
sample_metadata$mm_stage <- disease

#range(colSums(exprs(eset)))
# [1] 114848.0 117291.8

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# keep only healthy / myeloma samples
mask <- sample_metadata$disease %in% c('Healthy', 'Multiple Myeloma')

#table(mask)
# mask
# FALSE  TRUE 
#    31    25 

eset <- eset[, mask]
sample_metadata <- sample_metadata[mask, ]

# get gene symbols
gene_symbols <- fData(eset)$`Gene Symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

table(expr_dat$gene_symbol == '')
# 
# FALSE  TRUE 
# 20599  1000 

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

