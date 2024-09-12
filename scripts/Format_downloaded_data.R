# Format_downloaded_data.R

# This script formats and cleans clinical and expression data.
# - Creates "CLIN.txt" dimension 91 x 25
# - Creates 'expr_list.rds' including:
  # - expr_gene_tpm: dimension 61544 x 91
  # - expr_gene_counts: dimension 61544 x 91
  # - expr_isoform_tpm: dimension 246624 x 91
  # - expr_isoform_counts: dimension 246624 x 91

library(data.table)
library(readxl) 
library(stringr)
library(tximport)

# args <- commandArgs(trailingOnly = TRUE)
# work_dir <- args[1]
# annot_dir <- args[2]

# 1. CLIN_GIDE.txt

# Read clinical data
# clin <- read_excel(file.path(work_dir, '1-s2.0-S1535610819300376-mmc2.xlsx'), sheet='Table S1. PD-1 Patient')
clin_PD1 <- read_excel("files/1-s2.0-S1535610819300376-mmc2.xlsx")  # dim 56 x 15    # anti-PD-1 monotherapy: "Nivolumab" or "Pembrolizumab" treatment
clin_IPIPD1 <- read_excel("files/1-s2.0-S1535610819300376-mmc3.xlsx")  # dim 53 x 15    # combined anti-PSD-1 and anti-CTLA-4 : "Ipilimumab + Pembrolizumab" or "Ipilimumab + Nivolumab"
mapping <- read.csv("files/filereport_read_run_PRJEB23709_tsv.txt", sep = "\t")

# Table S1: PD-1, Patient 2, Pre-treatment = “PD1_2_PRE”.
# Table S2: IPIPD1, Patient 12, Pre-treatment = “ipiPD1_12_PRE”.

colnames(clin_PD1) <- clin_PD1[2, ] 
clin_PD1 <- clin_PD1[-c(1:2), ] 
colnames(clin_PD1) <- str_replace_all(colnames(clin_PD1), '\\W', '.')

colnames(clin_IPIPD1) <- clin_IPIPD1[2, ] 
clin_IPIPD1 <- clin_IPIPD1[-c(1:2), ] 
colnames(clin_IPIPD1) <- str_replace_all(colnames(clin_IPIPD1), '\\W', '.')

# Add "PD1_" prefix 
clin_PD1$Patient.no. <- paste0("PD1_", seq_len(nrow(clin_PD1)))
clin_IPIPD1$Patient.no. <- paste0("ipiPD1_", seq_len(nrow(clin_IPIPD1)))

# Bind the two clin data
clin <- rbind(clin_PD1, clin_IPIPD1)

# Split by semicolon and extractPatient.no. and treatment_timepoint 
mapping[, c("Patient.no.", "treatment_timepoint")] <- t(sapply(strsplit(mapping$submitted_ftp, "/"), function(x) {
  c(sub("(_PRE|_EDT).*", "", tail(x, 1)), sub(".*_(PRE|EDT).*", "\\1", tail(x, 1)))
}))

# Merge clinical data with mapping
clin_merged <- merge(clin, mapping, by = "Patient.no.")

# Rename run_accession to patient
colnames(clin_merged)[colnames(clin_merged) == "run_accession"] <- "patient"

# Save clinical data as CLIN.txt
write.table(clin_merged, "files/CLIN.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# write.table(clin, file.path(work_dir, "CLIN_GIDE.txt"), sep="\t", col.names=TRUE, row.names=FALSE)

# 2. expr_list.rds

# EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv
# source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
# load(file.path(annot_dir, "Gencode.v40.annotation.RData"))

# dir.create(file.path(work_dir, 'rnaseq'))
# zipfiles <- c('Gide_kallisto1.zip', 'Gide_kallisto2.zip', 'Gide_kallisto3.zip', 'Gide_kallisto4.zip', 'Gide_kallisto5.zip', 'Gide_kallisto6.zip')

source('https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/process_kallisto_output.R')
load("~/BHK lab/Annotation/Gencode.v40.annotation.RData")
work_dir <- "files/kallisto_v0.46.1_GRCh38.40"

# Creates EXPR_gene_tpm.tsv, EXPR_gene_counts.tsv, EXPR_tx_tpm.tsv, EXPR_tx_counts.tsv 
process_kallisto_output <- function(work_dir, tx2gene) {
  samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names = FALSE, recursive = FALSE)
  files <- file.path(work_dir, 'rnaseq', samples, "abundance.h5")
  names(files) <- samples
  
  expr_tx <- tximport(files, type = "kallisto", txOut = TRUE, ignoreAfterBar = TRUE)
  expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  
  expr_list <- list(
    expr_gene_tpm = log2(expr_gene$abundance + 0.001),
    expr_gene_counts = log2(expr_gene$counts + 1),
    expr_isoform_tpm = log2(expr_tx$abundance + 0.001),
    expr_isoform_counts = log2(expr_tx$counts + 1)
  )
  
  saveRDS(expr_list, file.path(work_dir, 'expr_list.rds'))
  unlink(file.path(work_dir, 'rnaseq'), recursive = TRUE)
}

# output: expr_list.rds
process_kallisto_output(work_dir, tx2gene)
