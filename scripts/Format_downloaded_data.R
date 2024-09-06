# Format_downloaded_data.R

# Formats and cleans clinical and expression data.
# Outputs: 
# - "CLIN.txt" (91 x ???)
# - 'expr_list.rds' including:
# - expr_gene_tpm: 61544 x 91
# - expr_gene_counts: 61544 x 43
# - expr_isoform_tpm: 246624 x 91
# - expr_isoform_counts: 246624 x 91

library(data.table)
library(readxl) 
library(stringr)
library(tximport)

# Clinical data
clin_PD1 <- read_excel("files/1-s2.0-S1535610819300376-mmc2.xlsx") # anti-PD-1 monotherapy (Treatmnest: "Nivolumab" or "Pembrolizumab")
clin_IPIPD1 <- read_excel("files/1-s2.0-S1535610819300376-mmc3.xlsx")  # ombined anti-PD-1 and anti-CTLA-4 immunotherapy (Treatmnest: "Ipilimumab + Pembrolizumab" or "Ipilimumab + Nivolumab")
mapping <- read.csv("files/filereport_read_run_PRJEB23709_tsv.txt", sep = "\t")

# Clean clinical data
colnames(clin_PD1) <- clin_PD1[2, ]
clin_PD1 <- clin_PD1[-c(1:2), ]
colnames(clin_IPIPD1) <- clin_IPIPD1[2, ]
clin_IPIPD1 <- clin_IPIPD1[-c(1:2), ]

# Add "PD1_" and "ipiPD1_" prefixes
clin_PD1$Patient.no. <- paste0("PD1_", seq_len(nrow(clin_PD1)))
clin_IPIPD1$Patient.no. <- paste0("ipiPD1_", seq_len(nrow(clin_IPIPD1)))

# Combine clinical data
clin <- rbind(clin_PD1, clin_IPIPD1)

# Extract Patient.no. from mapping
mapping$Patient.no. <- sapply(strsplit(mapping$submitted_ftp, "/"), function(x) {
  sub("(_PRE|_EDT).*", "", tail(x, n=1))
})

# Merge clinical and mapping data
clin_merged <- merge(clin, mapping, by = "Patient.no.")
colnames(clin_merged)[colnames(clin_merged) == "run_accession"] <- "patient"

# Save clinical data
write.table(clin_merged, "files/CLIN.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Load Gencode annotation
load("files/Gencode.v40.annotation.RData")

# Process expression data (kallisto output)
process_kallisto_output <- function(work_dir, tx2gene) {
  samples <- list.dirs(file.path(work_dir, 'rnaseq'), full.names = FALSE, recursive = FALSE)
  files <- file.path(work_dir, 'rnaseq', samples, "abundance.h5")
  names(files) <- samples
  
  # Import expression data
  expr_tx <- tximport(files, type = "kallisto", txOut = TRUE, ignoreAfterBar = TRUE)
  expr_gene <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  
  # Prepare expression data lists
  expr_list <- list(
    expr_gene_tpm = log2(expr_gene$abundance + 0.001),
    expr_gene_counts = log2(expr_gene$counts + 1),
    expr_isoform_tpm = log2(expr_tx$abundance + 0.001),
    expr_isoform_counts = log2(expr_tx$counts + 1)
  )
  
  return(expr_list)
}

# Define the working directory and process the expression data
work_dir <- "files/kallisto_v0.46.1_GRCh38.40"
expr_list <- process_kallisto_output(work_dir, tx2gene)

# Save the expression list as RDS
saveRDS(expr_list, file = "files/expr_list.rds")
