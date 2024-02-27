#!/bin/bash

#####
#
# MAIN CONFIGURATION
#

# Location of your FASTQ data
DATA_DIR='/mnt/data/bennettlab'
FASTQ_DIR="${DATA_DIR}/FASTQ"

# The locatiom of your reference files (Genome FASTA and indices etc.)
REF="/mnt/data/REFERENCES"

# The location of the PIPELINES
PIPELINE_DIR="${HOME}"
TRIMBA="${PIPELINE_DIR}/nf-core-trimba"

# The location of the nextflow exectuable
NEXTFLOW="nextflow"

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
GRCM38_GTF_DB="${REF}/nf-core/GRCm38/ensembl_gtf/default/Mus_musculus.GRCm38.102.gtf.db"
GRCM38_BOWTIE_IDX="${REF}/nf-core/GRCm38/fasta/default/GENOME/bowtie"
GRCM38_RNAS_BOWTIE_IDX="${REF}/nf-core/GRCm38/fasta/default/RNA/bowtie"
GRCM38_TRANSCRIPTS_BOWTIE_IDX="${REF}/nf-core/GRCm38/fasta/default/TRANSCRIPTOME/bowtie"
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
OUTDIR_TRIMBA="${DATA_DIR}/WF/trimba_demo"
OUTDIR_RNASEQ="${DATA_DIR}/WF/rnaseq_demo"

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

# Configuration for running on SLURM
cat > "${CONFIG_TRIMBA}" <<EOFCONFIGT
params {
  max_cpus = 100
  max_memory = 160.GB
}
process.executor = 'slurm'
process.queue = 'batch'
EOFCONFIGT

cat > "${CONFIG_RNASEQ}" <<EOFCONFIGR
params {
  max_cpus = 80
  max_memory = 128.GB
}
process.executor = 'slurm'
process.queue = 'bigmem'
EOFCONFIGR

#####
#
# TRIMBA_SAMPLESHEET - Edit this
#
# This is the main input file for the `trimba` pipeline
# If no riboseq data is available, only provide the header line.
cat > "${SAMPLESHEET_TRIMBA}" <<EOFSST
sample,fastq_1,fastq_2,strandedness
Riboseq_Control,"${FASTQ_DIR}/RIBO_DEMO/Ribo-Seq_of_MEF_Cell_Control.fastq.gz",,unstranded
Riboseq_Starvation,"${FASTQ_DIR}/RIBO_DEMO/Ribo-Seq_of_MEF_Cell_Amino_Acid_Starvation.fastq.gz",,unstranded
QTIseq_Control,"${FASTQ_DIR}/RIBO_DEMO/QTI-Seq_of_MEF_Cell_Control.fastq.gz",,unstranded
QTIseq_Starvation,"${FASTQ_DIR}/RIBO_DEMO/QTI-Seq_of_MEF_Cell_Amino_Acid_Starvation.fastq.gz",,unstranded
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
RNASeq_Control,"${FASTQ_DIR}/RIBO_DEMO/RNA-Seq_of_MEF_Cell_Control.fastq.gz",,auto
RNASeq_Starvation,"${FASTQ_DIR}/RIBO_DEMO/RNA-Seq_of_MEF_Cell_Amino_Acid_Starvation.fastq.gz",,auto
Riboseq_Control,"${OUTDIR_TRIMBA}/bowtie_rna/Riboseq_Control.unmapped.fastq.gz",,auto
Riboseq_Starvation,"${OUTDIR_TRIMBA}/bowtie_rna/Riboseq_Starvation.unmapped.fastq.gz",,auto
QTIseq_Control,"${OUTDIR_TRIMBA}/bowtie_rna/QTIseq_Control.unmapped.fastq.gz",,auto
QTIseq_Starvation,"${OUTDIR_TRIMBA}/bowtie_rna/QTIseq_Starvation.unmapped.fastq.gz",,auto
EOFSSR

#####
#
# GENES_OF_INTEREST - Edit this
#
# This contains the genes of interest for this `trimba` run
cat > "${GENES_OF_INTEREST}" <<EOFGOI
Brca1
Zfand2a
Tcea1
Gm7341
EOFGOI

#####
#
# MOTIFS_FILE - Edit this
#
# This contains the list of MOTIF files to be used in the known-motif part of the `trima` pipeline
# Note that you cannot mix DNA and RNA motif files, but for all RNA files there are DNA "equivalents"
cat > "${MOTIFS_FILE}" <<EOFMTFS
${MEME_DB}/CISBP-RNA/Mus_musculus.dna_encoded.meme
${MEME_DB}/CISBP-RNA/Homo_sapiens.dna_encoded.meme
${MEME_DB}/RNA/Ray2013_rbp_Mus_musculus.dna_encoded.meme
${MEME_DB}/RNA/Ray2013_rbp_Homo_sapiens.dna_encoded.meme
${MEME_DB}/MOUSE/HOCOMOCOv11_core_MOUSE_mono_meme_format.meme
${MEME_DB}/MOUSE/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme
${MEME_DB}/MOUSE/HOCOMOCOv10_MOUSE_mono_meme_format.meme
${MEME_DB}/MOUSE/chen2008.meme
${MEME_DB}/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
${MEME_DB}/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme
${MEME_DB}/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme
${MEME_DB}/HUMAN/HOCOMOCOv9.meme
${MEME_DB}/MIRBASE/22/Homo_sapiens_hsa.dna_encoded.meme
${MEME_DB}/MIRBASE/22/Mus_musculus_mmu.dna_encoded.meme
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
gtfdb:                          "$GRCM38_GTF_DB"
rna_index:                      "$GRCM38_RNAS_BOWTIE_IDX"
bowtie_index:                   "$GRCM38_BOWTIE_IDX"
transcripts_index:              "$GRCM38_TRANSCRIPTS_BOWTIE_IDX"

genes_file:                     "$GENES_OF_INTEREST"
motifs_file:                    "$MOTIFS_FILE"

skip_trimming:                  true
skip_riboseq:                   false
trimmer:                        "cutadapt"
extra_cutadapt_args:            '--trim-n --match-read-wildcards -u 16 -n 4 -a AGATCGGAAGAGCACACGTCTG -a AAAAAAAA -a GAACTCCAGTCAC -e 0.2 --nextseq-trim 20 -m 17 -M 34'
rnas_to_filter:                 "scRNA,3prime_overlapping_ncRNA,miRNA,snRNA,macro_lncRNA,sRNA,lincRNA,Mt_rRNA,scaRNA,snoRNA,rRNA,Mt_tRNA,bidirectional_promoter_lncRNA,misc_RNA"
nmotifs:                        1
extra_align_rna_args:           ''
extra_align_genome_args:        '-m 1'
extra_align_transcriptome_args: '-m 1'
extra_ribotish_qual_args:       ''
extra_ribotish_args:            ''

publish_dir_mode:               'symlink'
save_reference:                 false
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
salmon_index:          "${GRCM38_SALMON_INDEX}"
pseudo_aligner:        "salmon"

publish_dir_mode:      'symlink'
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