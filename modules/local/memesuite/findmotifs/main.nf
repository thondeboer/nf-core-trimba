process MEME_FINDMOTIFS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meme:5.5.4--pl5321hda358d9_0' :
        'ghcr.io/thondeboer/memesuite:5.5.5' }"
    // container "docker.io/memesuite/memesuite:5.5.5"
    // container "biocontainers/meme:5.5.4--pl5321hda358d9_0"
    // ghcr.io/thondeboer/memesuite:5.5.5

    input:
    tuple val(meta), path(fasta)
    val(nmotifs)

    output:
    tuple val(meta), path("*/*.eps"), emit: logoeps
    tuple val(meta), path("*/*.png"), emit: logopng, optional: true
    tuple val(meta), path("*/*.txt"), emit: resulttxt
    tuple val(meta), path("*/*.xml"), emit: resultxml
    tuple val(meta), path("*/*.html"), emit: htmlreport
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    meme \\
       -o ${prefix}_findmotifs \\
       -p $task.cpus \\
       -nmotifs $nmotifs \\
        $args \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meme: \$(meme -version)
    END_VERSIONS
    """
}
