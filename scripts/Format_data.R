# Format_Data Processing 
# Format_EXPR.R
# Goal: 
# - Creates "CLIN.txt" with dimensions 91 x 39
# - Creates cased_sequenced.csv: dimensions 91 x 4
# - Creates EXPR_gene_tpm.csv, EXPR_gene_counts.csv, EXPR_tx_tpm.csv, EXPR_tx_counts.csv 
#   - EXPR_gene_tpm.csv: dimensions 61,544 x 91
#   - EXPR_gene_counts.csv: dimensions 61,544 x 91
#   - EXPR_tx_tpm.csv: dimensions 246,624 x 91
#   - EXPR_tx_counts.csv: dimensions 246,624 x 91

# Libraries and Source Files
# Access necessary functions from ICB_common codes.
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/Get_Response.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/format_clin_data.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_tissue.R")
source("https://raw.githubusercontent.com/BHKLAB-Pachyderm/ICB_Common/main/code/annotate_drug.R")

library(data.table)
library(R.utils)
library(stringr)
library(tibble)

# args <- commandArgs(trailingOnly = TRUE)
# input_dir <- args[1]
# output_dir <- args[2]
# annot_dir <- args[3]

## 1. Get Clinical Data (CLIN.csv)
clin_orginal <- read.csv("files/CLIN.txt", stringsAsFactors=FALSE, sep="\t", dec=',') # 91 x 26

selected_cols <- c("patient" , "Age..Years." , "Sex" , "Treatment" , "Best.RECIST.response" , 
                   "Progression.Free.Survival..Days." , "Overall.Survival..Days." , "Progressed" , "Last.Followup.Status")

clin <- cbind(clin_orginal[, selected_cols], "Melanoma", NA, NA, NA, NA, NA, NA, NA, NA)

# Reorder columns
colnames(clin) <- c( "patient" , "age" , "sex" , "drug_type" , "recist" , "t.pfs" , "t.os" , "pfs" ,"os" ,
                     "primary" , "histo" , "response" , "stage" , "response.other.info" , "dna" ,"dna_info", "rna", "rna_info")

# Set rna and rna_info
clin$rna <- "rnaseq"
clin$rna_info <- "tpm"

# Set pfs and os, and convert OS days into months, ensuring two decimal points
clin$pfs <- as.numeric(as.character(ifelse(clin$pfs %in% "Yes", 1, ifelse(clin$pfs %in% "No", 0, NA))))
clin$t.pfs <- format(as.numeric(as.character(clin$t.pfs)) / 30)
clin$t.os <- format(as.numeric(as.character(clin$t.os)) / 30)

# Set os as binary representation
clin$os <- ifelse(clin$os %in% "Alive", 0, ifelse(clin$os %in% c("Dead", "Dead, melanoma"), 1, 0))

# Calculate the response using Get_Response function
clin$response <- Get_Response(data = clin)

# Reorder columns
clin <- clin[, c(
  "patient", "sex", "age", "primary", "histo", "stage", 
  "response.other.info", "recist", "response", "drug_type", "dna", "dna_info", "rna", "rna_info", "t.pfs", 
  "pfs", "t.os", "os"
)]

# Use the format_clin_data function for further formatting
clin <- format_clin_data(clin_orginal, "patient", selected_cols, clin)

# Read 'curation_tissue.csv' file
path <- "https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_tissue.csv"
annotation_tissue <- read.csv(path)

# Annotate 'clin' using the 'annotation_tissue' data; Set tissueid column (survival_unit and survival_type columns are added in this step)
clin <- annotate_tissue(clin=clin, study='Gide', annotation_tissue=annotation_tissue, check_histo=FALSE) 

# Set treatmentid after tissueid column, based on curation_drug.csv file
annotation_drug <- read.csv("https://raw.githubusercontent.com/BHKLAB-DataProcessing/ICB_Common/main/data/curation_drug.csv")
clin <- add_column(clin, treatmentid=annotate_drug('Gide', clin$drug_type, annotation_drug), .after='tissueid')

# Six categories: PD-1/PD-L1, CLA4, IO+combo, IO+chemo, Chemo+targeted, targeted, IO+targeted
clin$drug_type[clin$treatmentid == "Nivolumab" | clin$treatmentid == "Pembrolizumab"] <- 'PD-1/PD-L1'
clin$drug_type[clin$treatmentid == "Ipilimumab + Pembrolizumab" | clin$treatmentid == "Ipilimumab + Nivolumab"] <- 'IO+combo'

# Replace empty string values with NA
clin[clin == "-"] <- NA

# Remove any semicolons or unwanted characters in the columns
clin <- data.frame(lapply(clin, function(x) gsub(";", "", x)))

# Save the processed data as CLIN.csv file
write.table(clin, file="files/CLIN.csv", sep=";", quote=FALSE, col.names=TRUE, row.names=FALSE)

## 2. Get EXPR and case files
# Read expression data
expr_list <- readRDS("files/expr_list.rds")

# Patients should match between expression and clinical data
patients <- intersect(colnames(data.frame(expr_list$expr_gene_tpm)), clin$patient)
clin <- clin[clin$patient %in% patients, ]

# Create the case file
case <- data.frame(patient=patients, snv=0, cna=0, expr=1)
colnames(case) <- c("patient", "snv", "cna", "expr")

# Save the case file
write.table(case, file="files/cased_sequenced.csv", sep=";", quote=FALSE, col.names=TRUE, row.names=FALSE)

# Save the expression data
for (assay_name in names(expr_list)) {
  write.table(
    expr_list[[assay_name]],
    file=file.path("files/", paste0('EXPR_', str_replace(assay_name, 'expr_', ''), '.csv')),
    quote=FALSE, sep=";", col.names=TRUE, row.names=TRUE
  )
}

write.table( case , "files/cased_sequenced.csv" , quote=FALSE , sep=";" , col.names=TRUE , row.names=FALSE )
