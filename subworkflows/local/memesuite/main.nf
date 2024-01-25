include { MEME_FINDMOTIFS    } from '../../../modules/local/memesuite/findmotifs/main'
include { MEME_SEARCHMOTIFS  } from '../../../modules/local/memesuite/searchmotifs/main'

workflow MEMESUITE {
    take:
    fasta          // channel: [ val(meta), [ fastas ] ]
    motifs         // channel: [ motifs ]
    nmotifs        // int

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Find NOVEL motifs
    //
    ch_resulttxt = Channel.empty()
    MEME_FINDMOTIFS (
        fasta,
        nmotifs
    )
    ch_resulttxt = MEME_FINDMOTIFS.out.resulttxt
    ch_versions = ch_versions.mix(MEME_FINDMOTIFS.out.versions)

    //
    // MODULE: Find KNOWN motifs
    //
    MEME_SEARCHMOTIFS (
        fasta,
        motifs,
    )
    ch_sites = MEME_SEARCHMOTIFS.out.sites
    ch_versions = ch_versions.mix(MEME_SEARCHMOTIFS.out.versions)

    // TODO Add SEA output
    emit:
    resulttxt          = ch_resulttxt
    sites              = ch_sites
    versions           = ch_versions
}