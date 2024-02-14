#!/bin/bash

#####
#
# MAIN CONFIGURATION
#

# Location of your FASTQ data
DATA_DIR='/mnt/data/yasar'
FASTQ_DIR="${DATA_DIR}/FASTQ"

# The locatiom of your reference files (Genome FASTA and indices etc.)
REF="/mnt/data/REFERENCES"

# The location of the PIPELINES
PIPELINE_DIR="${HOME}/CODE/workflows"
TRIMBA="${PIPELINE_DIR}/nf-core-trimba"

# The location of the nextflow exectuable
NEXTFLOW="${HOME}/.local/bin/nextflow"

#####
#
# REFERENCE FILES CONFIGURATION
#

# Human data - REQUIRED - Only these THREE files are required
# See the documentation on how to download these files
GRCH38_FASTA="${REF}/nf-core/GRCh38/fasta/default/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GRCH38_GTF="${REF}/nf-core/GRCh38/ensembl_gtf/default/Homo_sapiens.GRCh38.110.gtf.gz"
MEME_DB="${REF}/memesuite/motif_databases"

# Mouse data - OPTIONAL (If not provided, will be created at run time but with save_reference can be saved for next time)
GRCH38_GTF_DB="${REF}/nf-core/GRCh38/ensembl_gtf/default/Homo_sapiens.GRCh38.110.gtf.db"
GRCH38_BOWTIE_IDX="${REF}/nf-core/GRCh38/fasta/default/GENOME/bowtie"
GRCH38_RNAS_BOWTIE_IDX="${REF}/nf-core/GRCh38/fasta/default/RNA/bowtie"
GRCH38_TRANSCRIPTS_BOWTIE_IDX="${REF}/nf-core/GRCh38/fasta/default/TRANSCRIPTOME/bowtie"
# For RNASEQ
GRCH38_STAR_INDEX="${REF}/nf-core/GRCh38/fasta/default/star"
GRCH38_SALMON_INDEX="${REF}/nf-core/GRCh38/fasta/default/salmon"

#*****************************
#* END OF MAIN CONFIGURATION *
#*****************************

#####
#
# OUTDIR CONFIGURATION - Edit this
#
# The OUTDIR determines where all the inputs and outputs are stored.
OUTDIR_TRIMBA="${DATA_DIR}/WF/trimba_firstrun"
OUTDIR_RNASEQ="${DATA_DIR}/WF/rnaseq_firstrun"

#####
#
# FILE SETUP - Do not edit
#
# This section just defines some more variables and should not be edited
mkdir -p $OUTDIR_TRIMBA
mkdir -p $OUTDIR_RNASEQ
SAMPLESHEET_TRIMBA="${OUTDIR_TRIMBA}/sample_sheet.csv"
SAMPLESHEET_RNASEQ="${OUTDIR_RNASEQ}/sample_sheet.csv"
GENES_OF_INTEREST="${OUTDIR_TRIMBA}/genes_of_interest.txt"
MOTIFS_FILE="${OUTDIR_TRIMBA}/motif_file.txt"
PARAMS_TRIMBA="${OUTDIR_TRIMBA}/nf-params.yml"
CONFIG_TRIMBA="${OUTDIR_TRIMBA}/nextflow.config"
PARAMS_RNASEQ="${OUTDIR_RNASEQ}/nf-params.yml"
CONFIG_RNASEQ="${OUTDIR_RNASEQ}/nextflow.config"

#####
#
# TRIMBA_SAMPLESHEET - Edit this
#
# This is the main input file for the `trimba` pipeline
# If no riboseq data is available, only provide the header line.
cat > "${SAMPLESHEET_TRIMBA}" <<EOFSST
sample,fastq_1,fastq_2,strandedness
EOFSST

#####
#
# RNASEQ_SAMPLESHEET - Edit this
#
# This contains the the input files for the `nf-core/rnaseq` pipeline.
# You should provide any orignal RNASEQ data and should also use the data in the `bowtie_rna` directory
# The latter will be created by the `trimba` pipeline, but you can anticipate what the name will be
cat > "${SAMPLESHEET_RNASEQ}" <<EOFSSR
sample,fastq_1,fastq_2,strandedness
EOFSSR

#####
#
# GENES_OF_INTEREST - Edit this
#
# This contains the genes of interest for this `trimba` run
cat > "${GENES_OF_INTEREST}" <<EOFGOI
Brca1
EOFGOI

#####
#
# MOTIFS_FILE - Edit this
#
# This contains the list of MOTIF files to be used in the known-motif part of the `trima` pipeline
# Note that you cannot mix DNA and RNA motif files, but for all RNA files there are DNA "equivalents"
cat > "${MOTIFS_FILE}" <<EOFMTFS
${MEME_DB}/CISBP-RNA/Mus_musculus.dna_encoded.meme
EOFMTFS

#####
#
# CONFIGURATION for running on SLURM
#
cat > "${CONFIG_TRIMBA}" <<EOFCONFIGR
params {
  max_cpus = 80
  max_memory = 128.GB
}
process.executor = 'slurm'
process.queue = 'bigmem'
process {
    withLabel:process_single {
        cpus   = 1
        memory = 4.GB
        time   = 4.h
    }
    withLabel:process_low {
        cpus   = 2
        memory = 8.GB
        time   = 8.h
    }
    withLabel:process_medium {
        cpus   = 4
        memory = 16.GB
        time   = 16.h
    }
    withLabel:process_high {
        cpus   = 8
        memory = 32.GB
        time   = 32.h
    }
}
EOFCONFIGR

cat > "${CONFIG_RNASEQ}" <<EOFCONFIGR
params {
  max_cpus = 80
  max_memory = 128.GB
}
process.executor = 'slurm'
process.queue = 'bigmem'
process {
    withLabel:process_single {
        cpus   = 1
        memory = 4.GB
        time   = 4.h
    }
    withLabel:process_low {
        cpus   = 2
        memory = 8.GB
        time   = 8.h
    }
    withLabel:process_medium {
        cpus   = 4
        memory = 16.GB
        time   = 16.h
    }
    withLabel:process_high {
        cpus   = 8
        memory = 31.GB
        time   = 31.h
    }
}
EOFCONFIGR

#####
#
# TRIMBA Pipeline parameters - Review and edit this (occasionally)
#
# These are the settings for the tools in the pipeline
# Review these, but you unlikley need to change these for different pipeline runs
cat > "${PARAMS_TRIMBA}" <<EOFPARAMST
input:                          "$SAMPLESHEET_TRIMBA"
outdir:                         "$OUTDIR_TRIMBA"
fasta:                          "$GRCH38_FASTA"
gtf:                            "$GRCH38_GTF"

genes_file:                     "$GENES_OF_INTEREST"
motifs_file:                    "$MOTIFS_FILE"

skip_trimming:                  true
skip_riboseq:                   true

save_reference:                 true
EOFPARAMST

#####
#
# RNASEQ Pipeline parameters - Review and edit this (occasionally)
#
# These are the settings for the tools in the pipeline
# Review these, but you unlikley need to change these for different pipeline runs
cat > "${PARAMS_RNASEQ}" <<EOFPARAMSR
input:                 "$SAMPLESHEET_RNASEQ"
outdir:                "$OUTDIR_RNASEQ"
fasta:                 "${GRCH38_FASTA}"
gtf:                   "${GRCH38_GTF}"
pseudo_aligner:        "salmon"
save_reference:        true
EOFPARAMSR

#####
#
# RNASEQ Extra options - Review and edit this (occasionally)
#
# The RNASEQ pipeline can be configured with MANY MORE options, and see the main nf-core/rnaseq website for options
# https://nf-co.re/rnaseq/3.14.0
# To run the most mimimal (and fastest) RNASEQ pipeline, you can skip a lot of the steps by providing these SKIP options
FAST="--skip_gtf_filter --skip_gtf_transcript_filter --skip_bbsplit --skip_umi_extract --skip_trimming \
      --skip_markduplicates --skip_bigwig --skip_stringtie --skip_fastqc --skip_preseq --skip_dupradar \
      --skip_qualimap --skip_rseqc --skip_biotype_qc --skip_deseq2_qc --skip_multiqc --skip_qc"
      
#************************************************************************************************************************
#* END OF CONFIGURATION                                                                                                 *
#************************************************************************************************************************

#####
#
# RUN TRIMBA pipeline
#
cd $OUTDIR_TRIMBA
$NEXTFLOW run $TRIMBA -profile docker -resume \
  -params-file "${PARAMS_TRIMBA}"
  
#####
#
# RUN RNASEQ pipeline
#
cd $OUTDIR_RNASEQ
$NEXTFLOW run nf-core/rnaseq -revision 3.14.0 -profile docker -resume \
  -params-file "${PARAMS_RNASEQ}" \
  $FAST