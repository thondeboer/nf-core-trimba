process COMBINETRIMBA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioframe:0.3.3--pyhdfd78af_0' :
        'quay.io/biocontainers/bioframe:0.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(transcript)
    tuple val(meta), path(ribotish)
    tuple val(meta), path(meme)
    tuple val(meta), path(sites)

    output:
    tuple val(meta), path("*.trimba_results.tsv")  , emit: trimbaresults, optional: true

    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    combinetrimba.py \\
        --transcript_file ${transcript} \\
        --ribotish_file ${ribotish} \\
        --meme_file ${meme} \\
        --sea_file ${sites} \\
        --prefix ${prefix}.trimba_results \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combinetrimba: \$(combinetrimba.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        combinetrimba: \$(combinetrimba.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}

