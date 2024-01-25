//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                          } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF                            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA               } from '../../../modules/nf-core/gunzip'
include { CUSTOM_GETCHROMSIZES                            } from '../../../modules/nf-core/custom/getchromsizes'
include { BOWTIE_BUILD as BOWTIE_RNA_BUILD                } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE_BUILD as BOWTIE_GENOME_BUILD             } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE_BUILD as BOWTIE_TRANSCRIPTOME_BUILD      } from '../../../modules/nf-core/bowtie/build'
include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA } from '../../../modules/nf-core/rsem/preparereference'


include { GETRNAS                             } from '../../../modules/local/getrnas'
include { GFF_CREATEDB                        } from '../../../modules/local/gff/createdb'

workflow PREPARE_GENOME {
    take:
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf
    bowtie_index         //      file: /path/to/genome.bowtie.idx
    gtfdb                //      file: /path/to/genome.gtf.db
    rna_index            //      file: /path/to/rna_bowtie_index
    rnas_to_filter       //      list: rnas to filter from genome fasta
    transcripts_fasta    //      file: /path/to/transcript.fasta
    transcripts_index    //      file: /path/to/transcript.fai
    save_reference       //      boolean: save reference genome files
    prepare_tool_indices //      list: tools to prepare indices for

    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(fasta))
    }

    //
    // Uncompress GTF annotation file if required
    //
    if (gtf) {
        if (gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(gtf))
        }
    } 

    //  Get GTF database
    if (!gtfdb) {
        ch_gtfdb    = GFF_CREATEDB ( [ [:], gtf ] ).db.map { it[1] }
        ch_versions = ch_versions.mix(GFF_CREATEDB.out.versions)
    }
    else {
        ch_gtfdb = Channel.value(file(gtfdb))
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    //
    // Create bowtie index of genome
    //
    if (!bowtie_index) {
        BOWTIE_GENOME_BUILD ( ch_fasta )
        ch_bowtie_index = BOWTIE_GENOME_BUILD.out.index
        ch_versions = ch_versions.mix( BOWTIE_GENOME_BUILD.out.versions )
    } else {
        ch_bowtie_index = Channel.value(file(bowtie_index))
    }
    //
    // Generate the bowtie indices from the selected RNAs or load from file
    //
    ch_rnas_index = Channel.empty()
    if ('RNAs' in prepare_tool_indices) {
        if (!rna_index) {
            //  Get RNAs
            ch_rna_fa = GETRNAS (
                ch_fasta.map { [ [:], it ] },
                ch_gtf,
                ch_gtfdb,
                rnas_to_filter
            ).fasta.map { it[1] }
            ch_versions = ch_versions.mix( GETRNAS.out.versions )
            // Create bowtie index of rnas

            BOWTIE_BUILD ( ch_rna_fa )
            ch_rnas_index = BOWTIE_BUILD.out.index
            ch_versions = ch_versions.mix( BOWTIE_BUILD.out.versions )
        } else {
            ch_rnas_index = Channel.value(file(rna_index))
        }
    } 

    //
    // Create bowtie index of transcriptome
    //
    ch_transcripts_fasta = Channel.empty()
    if (!transcripts_index) {
        //
        // Uncompress transcript fasta file / create if required
        //
        if (transcripts_fasta) {
            if (transcripts_fasta.endsWith('.gz')) {
                ch_transcripts_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], transcripts_fasta ] ).gunzip.map { it[1] }
                ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
            } else {
                ch_transcripts_fasta = Channel.value(file(transcripts_fasta))
            }
        } else {
            ch_transcripts_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_gtf ).transcripts_fasta
            ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
        }

        BOWTIE_TRANSCRIPTOME_BUILD ( ch_transcripts_fasta )
        ch_transcripts_index = BOWTIE_TRANSCRIPTOME_BUILD.out.index
        ch_versions = ch_versions.mix( BOWTIE_TRANSCRIPTOME_BUILD.out.versions )
    } else {
        ch_transcripts_index = Channel.value(file(transcripts_index))
    }

    emit:
    fasta             = ch_fasta                           // channel: path(genome.fasta)
    gtf               = ch_gtf                             // channel: path(genome.gtf)
    gtfdb             = ch_gtfdb                           // channel: path(genome.gtf.db)
    bowtie_index      = ch_bowtie_index                    // channel: path(bowtie/genome.bowtie.idx)
    rnas_index        = ch_rnas_index                      // channel: path(bowtie/rnas.index)
    fai               = ch_fai                             // channel: path(genome.fai)
    chrom_sizes       = ch_chrom_sizes                     // channel: path(genome.sizes)
    transcripts_fasta = ch_transcripts_fasta.ifEmpty(null) // channel: path(transcript.fasta)
    transcripts_index = ch_transcripts_index               // channel: path(bowtie/transcript.bowtie.idx)
    versions          = ch_versions.ifEmpty(null)          // channel: [ versions.yml ]
}
