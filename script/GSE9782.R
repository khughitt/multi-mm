#!/bin/env/Rscript
#
# Gene expression profiling and correlation with outcome in clinical trials of the proteasome inhibitor bortezomib
#
# Mulligan et al. (2007)
# n = 169 patients
#
# Description:
#
# The aims of this study were to assess the feasibility of prospective pharmacogenomics
# research in multicenter international clinical trials of bortezomib in multiple
# myeloma and to develop predictive classifiers of response and survival with
# bortezomib. Patients with relapsed myeloma enrolled in phase 2 and phase 3 clinical
# trials of bortezomib and consented to genomic analyses of pretreatment tumor samples.
# Bone marrow aspirates were subject to a negative-selection procedure to enrich for
# tumor cells, and these samples were used for gene expression profiling using DNA
# microarrays. Data quality and correlations with trial outcomes were assessed by
# multiple groups. Gene expression in this dataset was consistent with data published
# from a single-center study of newly diagnosed multiple myeloma. Response and survival
# classifiers were developed and shown to be significantly associated with outcome via
# testing on independent data. The survival classifier improved on the risk
# stratification provided by the International Staging System. Predictive models and
# biologic correlates of response show some specificity for bortezomib rather than
# dexamethasone. Informative gene expression data and genomic classifiers that predict
# clinical outcome can be derived from prospective clinical trials of new anticancer
# agents. Keywords: Gene expression profiling; correlation with outcome in clinical
# trials of the proteasome inhibitor bortezomib Overall design: Purified myeloma samples
# were collected prior to enrolment in clinical trials of bortezomib (PS-341). Samples
# were subject to replicate gene expression profiling using the Affymetrix 133A/B
# microarray. Data was normalized in MAS5.0 and the median of replicates is reported.
# Data was normalized to a Ttimmed mean of 15o and is NOT log transformed. Various
# patient parameters are reported as well as response, TTP and survival upon treatment
# with bortezomib or dexamethasone.
#
# Source: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA103715
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE9782'

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
# for GSE9782, data includes two separate esets with a small number of overlapping
# probes
esets <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)

# report data processing used
print(as.character(pData(esets[[1]])$data_processing[1]))
print(as.character(pData(esets[[2]])$data_processing[2]))

#
# Note: some of the metadata fields appear to be incorrectly encoded, e.g.:
#
# > pData(eset)[, 'characteristics_ch1.9']
# [1] PGx_Days_To_Progression = 87           PGx_Progression(0=No,1=Yes) = 1
# [3] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [5] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [7] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [9] PGx_CensorReason = No more data        PGx_Progression(0=No,1=Yes) = 1
#
# For the case of ch1.9, the column contains a mix of "pgx progression", 
# "days to progress", and "censor" fields. 
#

#
# ch1.7 - PGx_Response   (-> treatment_response)
# ch1.8 - PGx_Responder  (-> patient_subgroup)
#

# to get around this, we will manually detect and parse those fields..
sample_metadata <- pData(esets[[1]]) %>%
  mutate(
    study_code         = as.numeric(str_match(characteristics_ch1, '\\d+$')),
    treatment          = str_match(characteristics_ch1.1, '\\w+$'),
    gender             = str_match(characteristics_ch1.2, '\\w+$'),
    ethnicity          = str_match(characteristics_ch1.3, '\\w+$'),
    age                = as.numeric(str_match(characteristics_ch1.4, '\\d+$')),
    treatment_response = str_match(characteristics_ch1.7, '\\w+$'),
    patient_subgroup   = str_match(characteristics_ch1.8, '\\w+$')) %>%
  select(geo_accession, platform_id, study_code, treatment, gender,
          ethnicity, age, treatment_response, patient_subgroup) %>%
  add_column(geo_accession2 = pData(esets[[2]])$geo_accession, .after = 1)

# convert metadata from factor to character columns for parsing
str_mdat <- pData(esets[[1]])

for (cname in colnames(str_mdat)) {
  str_mdat[, cname] <- as.character(str_mdat[, cname])
}

pfs_event <- c()
pfs_time <- c()
pfs_event_reason <- c()
patient_died <- c()
os_time <- c()

# iterate over samples and find relevant information
for (sample_num in 1:nrow(pData(esets[[1]]))) {
  # PGx Progression
  ind <- which(grepl('PGx_Prog', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
    pfs_event <- c(pfs_event, 0)
  } else {
    pfs_event <- c(pfs_event, ifelse(endsWith(str_mdat[sample_num, ind], '1'), 1, 0))
  }

  # PGx days
  ind <- which(grepl('PGx_Days', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
    pfs_time <- c(pfs_time, NA)
  } else {
    pfs_time <- c(pfs_time, as.numeric(str_match(str_mdat[sample_num, ind], '\\d+$')))
  }

  # PGx Censor Reason
  ind <- which(grepl('PGx_CensorReason', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
      pfs_event_reason <- c(pfs_event_reason, NA)
  } else {
    pfs_event_reason <- c(pfs_event_reason, trimws(str_match(str_mdat[sample_num, ind], '[ \\w]+$')))
  }

  # Deceased
  ind <- which(grepl('Did_Patient_Die', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
    patient_died <- c(patient_died, 0)
  } else {
    patient_died <- c(patient_died, 1)
  }

  # Days Survived
  ind <- which(grepl('Days_Survived', str_mdat[sample_num, ]))
  os_time <- c(os_time, as.numeric(str_match(str_mdat[sample_num, ind], '\\d+$')))
}

# add to sample metadata
sample_metadata$pfs_time <- pfs_time
sample_metadata$pfs_event <- pfs_event
sample_metadata$pfs_event_reason <- pfs_event_reason
sample_metadata$os_time <- os_time
sample_metadata$patient_died <- patient_died

# add cell type and disease (same for all samples)
sample_metadata$disease = 'Multiple Myeloma'
sample_metadata$cell_type = NA # likely CD138+, but not 100% clear

# GSE9782 has *not* been adjusted for sample size..
#range(colSums(exprs(eset)))
# [1] 5772391 9324047

# for GSE7039, data includes two separate esets with a small number of overlapping
# probes
e1 <- exprs(esets[[1]])
e2 <- exprs(esets[[2]])

gene_symbols1 <- fData(esets[[1]])[, 'Gene symbol']
gene_symbols2 <- fData(esets[[2]])[, 'Gene symbol']

# adjust for size separately and then combine
e1 <- sweep(e1, 2, colSums(e1), '/') * 1E6
e2 <- sweep(e2, 2, colSums(e2), '/') * 1E6

mask1 <- !startsWith(rownames(e1), 'AFFX-')
mask2 <- !startsWith(rownames(e2), 'AFFX-')

e1 <- e1[mask1, ]
e2 <- e2[mask2, ]

gene_symbols1 <- gene_symbols1[mask1]
gene_symbols2 <- gene_symbols2[mask2]

# small number of shared probes exist...
#length(intersect(rownames(e1), rownames(e2)))
# [1] 168

#head(intersect(rownames(e1), rownames(e2)))
# [1] "200000_s_at" "200001_at"   "200002_at"   "200003_s_at" "200004_at"   "200005_at"

# use average values for those probes..
shared_probes <- intersect(rownames(e1), rownames(e2))

e1[shared_probes, ] <- (e1[shared_probes, ] + e2[shared_probes, ]) / 2

mask <- !rownames(e2) %in% shared_probes
e2 <- e2[mask, ]
gene_symbols2 <- gene_symbols2[mask]

expr_dat <- rbind(e1, e2)
gene_symbols <- c(gene_symbols1, gene_symbols2)

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

