# Need to be tested
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"], 
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)
prefix = config["prefix"]
filename = config["filename"]

# Rule to generate the MultiAssayExp object using the get_MultiAssayExp.R script
rule get_MultiAssayExp:
    input:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv"),
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + filename)
    resources:
        mem_mb=3000,
        disk_mb=3000
    shell:
        """
        Rscript /mnt/data/get_MultiAssayExp.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation \
        {prefix}{filename}
        """

# Rule to format data using Format_data.R
rule format_data:
    input:
        S3.remote(prefix + "download/expr_list.rds"),
        S3.remote(prefix + "download/mel_gide19_survival_data.csv"),
        S3.remote(prefix + "download/CLIN_GIDE.txt"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    output:
        S3.remote(prefix + "processed/cased_sequenced.csv"),
        S3.remote(prefix + "processed/CLIN.csv"),
        S3.remote(prefix + "processed/EXPR_gene_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_gene_counts.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_tpm.csv"),
        S3.remote(prefix + "processed/EXPR_isoform_counts.csv")
    shell:
        """
        Rscript scripts/Format_data.R \
        {prefix}download \
        {prefix}processed \
        {prefix}annotation
        """

# Rule to format downloaded data using Format_downloaded_data.R
rule format_downloaded_data:
    input:
        S3.remote(prefix + "download/Gide_kallisto1.zip"),
        S3.remote(prefix + "download/Gide_kallisto2.zip"),
        S3.remote(prefix + "download/Gide_kallisto3.zip"),
        S3.remote(prefix + "download/Gide_kallisto4.zip"),
        S3.remote(prefix + "download/Gide_kallisto5.zip"),
        S3.remote(prefix + "download/Gide_kallisto6.zip"),
        S3.remote(prefix + "download/1-s2.0-S1535610819300376-mmc2.xlsx"),
        S3.remote(prefix + "download/1-s2.0-S1535610819300376-mmc3.xlsx"),  # Added this file
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData")
    output:
        S3.remote(prefix + "download/CLIN_GIDE.txt"),
        S3.remote(prefix + "download/expr_list.rds")
    shell:
        """
        Rscript scripts/format_downloaded_data.R \
        {prefix}download \
        {prefix}annotation \
        """

# Rule to download annotation files
rule download_annotation:
    output:
        S3.remote(prefix + "annotation/Gencode.v40.annotation.RData"),
        S3.remote(prefix + "annotation/curation_drug.csv"),
        S3.remote(prefix + "annotation/curation_tissue.csv")
    shell:
        """
        wget https://github.com/BHKLAB-Pachyderm/Annotations/blob/master/Gencode.v40.annotation.RData?raw=true -O {prefix}annotation/Gencode.v40.annotation.RData 
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_drug.csv -O {prefix}annotation/curation_drug.csv
        wget https://github.com/BHKLAB-Pachyderm/ICB_Common/raw/main/data/curation_tissue.csv -O {prefix}annotation/curation_tissue.csv
        """

# Rule to download data files
rule download_data:
    output:
        S3.remote(prefix + "download/Gide_kallisto1.zip"),
        S3.remote(prefix + "download/Gide_kallisto2.zip"),
        S3.remote(prefix + "download/Gide_kallisto3.zip"),
        S3.remote(prefix + "download/Gide_kallisto4.zip"),
        S3.remote(prefix + "download/Gide_kallisto5.zip"),
        S3.remote(prefix + "download/Gide_kallisto6.zip"),
        S3.remote(prefix + "download/mel_gide19_survival_data.csv"),
        S3.remote(prefix + "download/1-s2.0-S1535610819300376-mmc2.xlsx"),
        S3.remote(prefix + "download/1-s2.0-S1535610819300376-mmc3.xlsx"),  # Added this file
        S3.remote(prefix + "download/filereport_read_run_PRJEB23709_tsv.txt")  # Added this file as well
    shell:
        """
        wget -O {prefix}download/Gide_kallisto1.zip "https://zenodo.org/record/6968597/files/Gide_kallisto1.zip?download=1"
        wget -O {prefix}download/Gide_kallisto2.zip "https://zenodo.org/record/6968597/files/Gide_kallisto2.zip?download=1"
        wget -O {prefix}download/Gide_kallisto3.zip "https://zenodo.org/record/6968597/files/Gide_kallisto3.zip?download=1"
        wget -O {prefix}download/Gide_kallisto4.zip "https://zenodo.org/record/6968597/files/Gide_kallisto4.zip?download=1"
        wget -O {prefix}download/Gide_kallisto5.zip "https://zenodo.org/record/6968597/files/Gide_kallisto5.zip?download=1"
        wget -O {prefix}download/Gide_kallisto6.zip "https://zenodo.org/record/6968597/files/Gide_kallisto6.zip?download=1"
        wget https://ars.els-cdn.com/content/image/1-s2.0-S1535610819300376-mmc2.xlsx -O {prefix}download/1-s2.0-S1535610819300376-mmc2.xlsx
        wget https://ars.els-cdn.com/content/image/1-s2.0-S1535610819300376-mmc3.xlsx -O {prefix}download/1-s2.0-S1535610819300376-mmc3.xlsx
        wget "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB23709&result=read_run&format=tsv&download=true" -O {prefix}download/filereport_read_run_PRJEB23709_tsv.txt
        """
