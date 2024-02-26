/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet} from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowTrimba.initialise(params, log)

def prepareToolIndices  = []
// Check rRNA databases for remove rRNA option
if (!params.skip_remove_noncoding_rna) {
    if (!params.fasta && !params.gtf  && !params.rna_index) { exit 1, "ERROR: You must provide a FASTA and GTF file to remove non-coding RNA, or provide a pre-made index, or use the --skip_remove_noncoding_rna flag to skip this step." }
    prepareToolIndices << 'RNAs'
}

ch_motifs_file = file(params.motifs_file)
if (ch_motifs_file.isEmpty()) {exit 1, "File provided with --motifs_file is empty: ${ch_motifs_file.getName()}!"}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
def dummy_file1                       = "$baseDir/assets/dummy_file.txt"
def dummy_file2                       = "$baseDir/assets/dummy_file2.txt"
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_GENOME   } from '../subworkflows/local/prepare_genome'
include { RIBOTISH         } from '../subworkflows/local/ribotish/main'
include { MEMESUITE        } from '../subworkflows/local/memesuite/main'

//
// MODULE: Installed directly from local
//
include { GETTRANSCRIPTS   } from '../modules/local/gettranscripts/main'
include { COMBINETRIMBA    } from '../modules/local/combinetrimba/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                                  } from '../modules/nf-core/cat/fastq'
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { FASTQC as POSTTRIM_FASTQC                  } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                                   } from '../modules/nf-core/cutadapt/main'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_RNA           } from '../modules/nf-core/bowtie/align'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_GENOME        } from '../modules/local/bowtie/align'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_TRANSCRIPTOME } from '../modules/local/bowtie/align'
include { BEDTOOLS_GENOMECOV                         } from '../modules/nf-core/bedtools/genomecov/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []
def pass_mapped_reads  = [:]
def pass_trimmed_reads = [:]
def pass_strand_check  = [:]

workflow TRIMBA {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME (
        params.fasta,
        params.gtf,
        params.bowtie_index,
        params.gtfdb,
        params.rna_index,
        params.rnas_to_filter,
        params.transcripts_fasta,
        params.transcripts_index,
        params.save_reference,
        prepareToolIndices,
    )
    // Get the bowtie1 index (Bowtie 1 is better/faster for the short reads we have)
    ch_bowtie_rnas_index = PREPARE_GENOME.out.rnas_index
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // Check if we even need to run riboSEQ
    //
    // Channel
    //     .fromPath(dummy_file)
    //     .into { dummyFileChan1; dummyFileChan2 }
    ch_ribotish_results = [ [:], dummy_file1 ]
    ch_genomecov        = [ [:], dummy_file2 ]
    if (!params.skip_riboseq) {

        //
        // Create input channel from input file provided through params.input
        //
        Channel
            .fromSamplesheet("input")
            .map {
                meta, fastq_1, fastq_2 ->
                    if (!fastq_2) {
                        return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                    } else {
                        return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                    }
            }
            .groupTuple()
            .map {
                WorkflowTrimba.validateInput(it)
            }
            .branch {
                meta, fastqs ->
                    single  : fastqs.size() == 1
                        return [ meta, fastqs.flatten() ]
                    multiple: fastqs.size() > 1
                        return [ meta, fastqs.flatten() ]
            }
            .set { ch_fastq }

        //
        // MODULE: Concatenate FastQ files from same sample if required
        //
        CAT_FASTQ (
            ch_fastq.multiple
        )
        .reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))


        //
        // Many options for trimming, so we have a few subworkflows
        //
        ch_filtered_reads      = ch_cat_fastq
        ch_fastqc_raw_multiqc  = Channel.empty()
        ch_fastqc_trim_multiqc = Channel.empty()
        ch_trim_log_multiqc    = Channel.empty()
        ch_trim_read_count     = Channel.empty()

        //
        // MODULE: CUTADAPT
        //
        // We mimic the subworkflow here for cutadapt with before and after trimming FASTQC
        if (params.trimmer == 'cutadapt' && !params.skip_trimming) {
            if (!params.skip_qc && !params.skip_fastqc) {
                FASTQC (
                    ch_cat_fastq
                )
                ch_fastqc_raw_multiqc = FASTQC.out.zip
                ch_versions = ch_versions.mix(FASTQC.out.versions.first())
            }
            CUTADAPT ( ch_cat_fastq )
            ch_filtered_reads      = CUTADAPT.out.reads
            ch_versions = ch_versions.mix(CUTADAPT.out.versions)
            if (!params.skip_qc && !params.skip_fastqc) {
                POSTTRIM_FASTQC (
                    ch_filtered_reads
                )
                ch_fastqc_trim_multiqc = POSTTRIM_FASTQC.out.zip
                ch_versions = ch_versions.mix(POSTTRIM_FASTQC.out.versions.first())
            }
            // TODO get trim_read_count from cutadapt
        }
        //
        // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
        //
        if (params.trimmer == 'trimgalore') {
            FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
                ch_cat_fastq,
                params.skip_fastqc || params.skip_qc,
                params.with_umi,
                params.skip_umi_extract,
                params.skip_trimming,
                params.umi_discard_read,
                params.min_trimmed_reads
            )
            ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
            ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip
            ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip
            ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log
            ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_read_count
            ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)
        }

        //
        // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
        //
        if (params.trimmer == 'fastp') {
            FASTQ_FASTQC_UMITOOLS_FASTP (
                ch_cat_fastq,
                params.skip_fastqc || params.skip_qc,
                params.with_umi,
                params.skip_umi_extract,
                params.umi_discard_read,
                params.skip_trimming,
                [],
                params.save_trimmed,
                params.save_trimmed,
                params.min_trimmed_reads
            )
            ch_filtered_reads      = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
            ch_fastqc_raw_multiqc  = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_raw_zip
            ch_fastqc_trim_multiqc = FASTQ_FASTQC_UMITOOLS_FASTP.out.fastqc_trim_zip
            ch_trim_log_multiqc    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_json
            ch_trim_read_count     = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_read_count
            ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
        }

        //
        // Get list of samples that failed trimming threshold for MultiQC report
        //
        ch_trim_read_count
            .map {
                meta, num_reads ->
                    pass_trimmed_reads[meta.id] = true
                    if (num_reads <= params.min_trimmed_reads.toFloat()) {
                        pass_trimmed_reads[meta.id] = false
                        return [ "$meta.id\t$num_reads" ]
                    }
            }
            .collect()
            .map {
                tsv_data ->
                    def header = ["Sample", "Reads after trimming"]
                    WorkflowTrimba.multiqcTsvFromList(tsv_data, header)
            }
            .set { ch_fail_trimming_multiqc }



        //
        // MODULE: Filter with bowtie1 to remove rRNA
        //
        if (!params.skip_remove_noncoding_rna) {
            BOWTIE_ALIGN_RNA ( 
                ch_filtered_reads,
                ch_bowtie_rnas_index
            )
            ch_versions = ch_versions.mix(BOWTIE_ALIGN_RNA.out.versions)
            ch_norna_reads = BOWTIE_ALIGN_RNA.out.fastq
        }

        //
        // MODULE: Align reads to genome
        //
        BOWTIE_ALIGN_GENOME (
            ch_norna_reads,
            PREPARE_GENOME.out.bowtie_index
        )
        ch_versions = ch_versions.mix(BOWTIE_ALIGN_GENOME.out.versions)

        //
        // MODULE: Align reads to transcriptome
        //
        BOWTIE_ALIGN_TRANSCRIPTOME (
            ch_norna_reads,
            PREPARE_GENOME.out.transcripts_index
        )
        ch_versions = ch_versions.mix(BOWTIE_ALIGN_TRANSCRIPTOME.out.versions)

        ch_bam_for_genomecov = BOWTIE_ALIGN_TRANSCRIPTOME.out.bam.map { meta, bam ->
            return [meta, bam, 1]
        }
        BEDTOOLS_GENOMECOV (
            ch_bam_for_genomecov,
            PREPARE_GENOME.out.chrom_sizes,
            'bedgraph'       
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)    //
        // SUBWORKFLOW: RiboTish Quality and Predict
        //
        RIBOTISH (
            BOWTIE_ALIGN_GENOME.out.bam,
            BOWTIE_ALIGN_GENOME.out.bai,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.fasta,
        )
        ch_versions = ch_versions.mix(RIBOTISH.out.versions)

        ch_ribotish_results = RIBOTISH.out.results
        ch_genomecov        = BEDTOOLS_GENOMECOV.out.genomecov
    } else {
        // We are skipping riboSEQ, so we need to create dummy channels
        def meta = [id: 'gene_list']
        ch_ribotish_results = [ meta, dummy_file1 ]
        ch_genomecov        = [ meta, dummy_file2 ]        
    }

    //
    // MODULE: Motif finding
    //
    def genes_file = { params.genes_file ?: dummy_file }
    GETTRANSCRIPTS (
        ch_ribotish_results,
        ch_genomecov,
        PREPARE_GENOME.out.gtfdb,
        PREPARE_GENOME.out.fasta,
        genes_file,
    )

    //
    // MODULE: Find known and novel motifs
    //
    ch_meme_results = [ [:], dummy_file1 ]
    ch_sites = [ [:], dummy_file2 ]
    if (!params.skip_meme) {
        ch_motifs= Channel.from(ch_motifs_file.readLines()).map { row -> file(row, checkIfExists: true) }.collect()
        MEMESUITE (
            GETTRANSCRIPTS.out.fasta,
            ch_motifs,
            params.nmotifs
        )
        ch_versions = ch_versions.mix(MEMESUITE.out.versions)
        ch_meme_results = MEMESUITE.out.resulttxt
        ch_sites = MEMESUITE.out.sites
    }
    //
    // MODULE: Combine all the results
    //
    COMBINETRIMBA (
        GETTRANSCRIPTS.out.fasta,
        GETTRANSCRIPTS.out.ribotish,
        ch_meme_results,
        ch_sites,
    )
    ch_versions = ch_versions.mix(COMBINETRIMBA.out.versions)


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    
    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc & !params.skip_riboseq) {
        workflow_summary    = WorkflowTrimba.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowTrimba.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_raw_multiqc.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_trim_multiqc.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
