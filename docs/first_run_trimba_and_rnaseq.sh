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

# Mouse data - REQUIRED - Only these THREE files are required
# See the documentation on how to download these files
GRCM38_FASTA="${REF}/nf-core/GRCm38/fasta/default/Mus_musculus.GRCm38.dna.primary_assembly.fa"
GRCM38_GTF="${REF}/nf-core/GRCm38/ensembl_gtf/default/Mus_musculus.GRCm38.102.gtf.gz"
MEME_DB="${REF}/memesuite/motif_databases"

# Mouse data - OPTIONAL (If not provided, will be created at run time but with save_reference can be saved for next time)
GRCM38_GTF_DB="${REF}/nf-core/GRCm38/ensembl_gtf/default/gtf.db"
GRCM38_BOWTIE_IDX="${REF}/nf-core/GRCm38/fasta/default/bowtie_index"
GRCM38_RNAS_BOWTIE_IDX="${REF}/nf-core/GRCm38/fasta/default/Mus_musculus.GRCm38.dna.primary_assembly.RNAtranscripts_bowtie"
GRCM38_TRANSCRIPTS="${REF}/nf-core/GRCm38/fasta/default/Mus_musculus.GRCm38.dna.primary_assembly.ALLtranscripts.fa"
GRCM38_TRANSCRIPTS_BOWTIE_IDX="${REF}/nf-core/GRCm38/fasta/default/Mus_musculus.GRCm38.dna.primary_assembly.ALLtranscripts_bowtie"
# For RNASEQ
GRCM38_STAR_INDEX="${REF}/nf-core/GRCm38/fasta/default/star"
GRCM38_SALMON_INDEX="${REF}/nf-core/GRCm38/fasta/default/salmon"

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
# TRIMBA Pipeline parameters - Review and edit this (occasionally)
#
# These are the settings for the tools in the pipeline
# Review these, but you unlikley need to change these for different pipeline runs
cat > "${PARAMS_TRIMBA}" <<EOFPARAMST
input:                          "$SAMPLESHEET_TRIMBA"
outdir:                         "$OUTDIR_TRIMBA"
fasta:                          "$GRCM38_FASTA"
gtf:                            "$GRCM38_GTF"

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
fasta:                 "${GRCM38_FASTA}"
gtf:                   "${GRCM38_GTF}"
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
FAST="--skip_gtf_filter --skip_gtf_transcript_filter --skip_bbsplit --skip_umi_extract --skip_trimming --skip_alignment \
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