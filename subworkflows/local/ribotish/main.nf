include { RIBOTISH_QUALITY } from '../../../modules/local/ribotish/quality/main'
include { RIBOTISH_PREDICT } from '../../../modules/local/ribotish/predict/main'

workflow RIBOTISH {
    take:
    bam          // channel: [ val(meta), [ bams ] ]
    bai          // channel: [ val(meta), [ bais ] ]
    gtf          // path gtf
    fasta        // path fasta

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: RiboTish Quality
    //
    RIBOTISH_QUALITY (
        bam,
        bai,
        gtf,
    )
    ch_versions = ch_versions.mix(RIBOTISH_QUALITY.out.versions)

    //
    // MODULE: RiboTish
    //
    RIBOTISH_PREDICT (
        bam,
        bai,
        RIBOTISH_QUALITY.out.psites,
        gtf,
        fasta,
    )
    ch_versions = ch_versions.mix(RIBOTISH_PREDICT.out.versions)

    emit:
    results            = RIBOTISH_PREDICT.out.results
    allresults         = RIBOTISH_PREDICT.out.allresults
    versions           = ch_versions


}