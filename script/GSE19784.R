#!/bin/env/Rscript
#
# Gene expression profiling for molecular classification of multiple myeloma in newly
# diagnosed patients
# Broyl et al. (2010)
#
# Clinical trial: http://www.hovon.nl/studies/studies-per-ziektebeeld/mm.html?action=showstudie&studie_id=5&categorie_id=3
# 
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE19784'

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
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# report data processing used
print(as.character(pData(eset)$data_processing[1]))

# exclude control sequences present in some datasets (GSE19784)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# size factor scaling has already been performed
#range(colSums(exprs(eset))) 
# [1] 212082.1 239655.9

# in order to normalize downstream comparisons across datasets, however, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

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

#
# fields of interest:
#
# cluster:ch1                                                                                                                                                                                    
# CD-1       CD-2        CTA         HY         MF         MS    Myeloid       NFκB       None                                                                                             
# 13         34         22         80         32         33         40         39          9                                                                                             
# PR       PRL3 SOCS3/PRL3                                                                                                                                                               
# 15          2          9                                                                                                                                                               
#iss:ch1                                                                                                                                                                                        
#  I  II III  nd                                                                                                                                                                                
#122  88  83  27        
#

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(sample_id = geo_accession, platform_id, iss_stage = `iss:ch1`,
         patient_subgroup = `cluster:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'BM-CD138+'

# load additional metadata provided in supplemental table S11
# Note; there is not sufficient information provided to link patients in table S11
# to GSM sample identifiers..
#table_s11 <- read_csv(file.path(base_dir, 'metadata', 'broyl2010_supp_table_s11.csv'))

survival_mdata <- read_csv(file.path(base_dir, 'metadata', 'kuiper2012_supp_patient_survival.csv'))

survival_mdata <- survival_mdata %>%
  rename(sample_id = Patient) %>%
  filter(sample_id %in% colnames(eset))

colnames(survival_mdata) <- c('sample_id', 'os_months', 'os_event', 'pfs_months', 'pfs_event')

# exclude samples without metadata
mask <- sample_metadata$sample_id %in% survival_mdata$sample_id
 
#all(colnames(eset) == sample_metadata$sample_id)
# [1] TRUE

eset <- eset[, mask]
sample_metadata <- sample_metadata[mask, ]

# combine metadata
sample_metadata <- sample_metadata %>%
  inner_join(survival_mdata, by = 'sample_id')

# load expression data and add gene symbol column
expr_dat <- as.data.frame(exprs(eset))

expr_dat <- expr_dat %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = fData(eset)$`Gene symbol`, .after = 1)

# store cleaned expression data and metadata
write_csv(expr_dat, path = file.path(clean_data_dir, sprintf('%s_1_expr.csv', accession)))
write_csv(sample_metadata, 
          file.path(clean_data_dir, sprintf('%s_1_sample_metadata.csv', accession)))

sessionInfo()

