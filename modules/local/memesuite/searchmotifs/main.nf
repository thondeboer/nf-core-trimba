process MEME_SEARCHMOTIFS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.4--pl5321hda358d9_0' :
        'ghcr.io/thondeboer/memesuite:5.5.5' }"
    // container "docker.io/memesuite/memesuite:5.5.5"
    // container "biocontainers/meme:5.5.4--pl5321hda358d9_0"
    // quay.io/biocontainers/meme:5.5.5--pl5321hda358d9_0
    // ghcr.io/thondeboer/memesuite:5.5.5

    input:
    tuple val(meta), path(fasta)
    path motifs

    output:
    tuple val(meta), path("*/sea.html")       , emit: htmlreport, optional: true
    tuple val(meta), path("*/sea.tsv")        , emit: resulttsv , optional: true
    tuple val(meta), path("*/sequences.tsv")  , emit: sequences , optional: true
    tuple val(meta), path("*/sites.tsv")      , emit: sites     , optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sea \\
       ${args} \\
       --o ${prefix} \\
       --p ${fasta} \\
       ${'--m '+motifs.join(' --m ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meme: \$(sea -version)
    END_VERSIONS
    """
}
