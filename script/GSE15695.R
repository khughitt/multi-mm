#!/bin/env/Rscript
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE15695'

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
esets <- getGEO(accession, destdir = raw_data_dir)

# GSE15695 includes a mix of SNP array and microarray data;
# the GPL with the microarray data is GPL570

# lapply(esets, dim)                                                                         
#$`GSE15695-GPL3718_series_matrix.txt.gz`                                                     
#Features  Samples                                                                            
#  262314      168                                                                            
                                                                                             
#$`GSE15695-GPL3720_series_matrix.txt.gz`                                                     
#Features  Samples                                                                            
#  238354      168                                                                            
                                                                                             
#$`GSE15695-GPL570_series_matrix.txt.gz`                                                      
#Features  Samples                                                                            
#   54675      247     

eset <- esets[[3]]

# report data processing used
print(pData(eset)$data_processing[1])

# GSE15695 has not been size-factor scaled
#range(colSums(exprs(eset)))
# [1] 13401102 26073325

# in order to normalize downstream comparisons across datasets, however, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets (GSE6477)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'BM-CD138+'

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = fData(eset)$`Gene Symbol`, .after = 1)

# store cleaned expression data and metadata
expr_outfile <- file.path(clean_data_dir, sprintf('%s_gene_expr.csv', accession))
mdat_outfile <- file.path(clean_data_dir, sprintf('%s_sample_metadata.csv', accession))

write_csv(expr_dat, expr_outfile)
write_csv(sample_metadata, mdat_outfile)

sessionInfo()

