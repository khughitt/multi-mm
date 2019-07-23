#!/bin/env/Rscript
#
#	Molecular prognosis in multiple myeloma (2008)
# n = 250 patients
#
# Note: Survival data provided by corresponding author (Stéphane Minvielle) via email 
# on June 6, 2019.
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE7039'

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
# for GSE7039, data includes two separate esets with a small number of overlapping
# probes; Also note that the data corresponds to two platforms, only one of which
# has a corresponding .annot.gz file.

# The GEO project includes two replicates ("_A" and "_B") for each patient.
esets <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)

# report data processing used
print(as.character(pData(esets[[1]])$data_processing[1]))

# load additional survival metadata
survival_dat <- read_csv(file.path(base_dir, 'metadata', 'MM survival time GSE7039.csv'))

# combine samples from separate ExpressionSets
sample_metadata <- pData(esets[[1]]) %>%
  mutate(patient_id = sub('_A', '', title)) %>%
  select(geo_accession, platform_id, patient_id) %>%
  inner_join(survival_dat, by = 'patient_id') %>%
  add_column(geo_accession2 = pData(esets[[2]])$geo_accession, .after = 1) %>%
  mutate(patient_died = ifelse(deceased == 'yes', 1, 0)) %>%
  rename(os_time = follow_up_days) %>%
  select(-deceased)

#all(sample_metadata$patient_id == sub('_B', '', pData(esets[[2]])$title))
# [1] TRUE

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = 'CD138+'
sample_metadata$tissue = 'Bone Marrow'
sample_metadata$sample_type = 'Patient'

# adjust for size separately and then combine
e1 <- exprs(esets[[1]])
e2 <- exprs(esets[[2]])

e1 <- sweep(e1, 2, colSums(e1), '/') * 1E6
e2 <- sweep(e2, 2, colSums(e2), '/') * 1E6

mask1 <- !startsWith(rownames(e1), 'AFFX-')
mask2 <- !startsWith(rownames(e2), 'AFFX-')

e1 <- e1[mask1, ]
e2 <- e2[mask2, ]

# no shared probe ids
# length(intersect(rownames(e1), rownames(e2)))
# [1] 0
expr_dat <- rbind(e1, e2)

# get gene symbols
gene_symbols <- c(fData(esets[[1]])[, 'Gene symbol'][mask1],
                  fData(esets[[2]])[, 'Gene Symbol'][mask2])

# drop ambiguous / non-gene fields
# *Multi Hs   *Genomic sequence *Repeats containing   *Seq not verified                ESTs                                                                                            
#             644                 577                 567                 246                 162                                                                                            
mask <- !startsWith(gene_symbols, '*')

#table(mask)
# mask
# FALSE  TRUE 
#  2156 10187 

expr_dat <- expr_dat[mask, ]
gene_symbols <- gene_symbols[mask]

# get expression data and add gene symbol column
expr_dat <- expr_dat %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(gene_symbol = gene_symbols, .after = 1)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.csv', accession)
mdat_outfile <- sprintf('%s_sample_metadata.csv', accession)

# store cleaned expression data and metadata
write_csv(expr_dat, file.path(clean_data_dir, expr_outfile))
write_csv(sample_metadata, file.path(clean_data_dir, mdat_outfile))

sessionInfo()

