process GETRNAS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffutils:0.11.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/gffutils:0.12--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(fasta)
    path gtf
    path dbfile
    val rnas

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val(meta), path("*.db"), emit: gtfdb, optional: true

    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    get_rnas.py \\
        --fasta $fasta \\
        --gtf $gtf \\
        --rnas_to_filter $rnas \\
        --db_file $dbfile \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getrnas: \$(get_rnas.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getrnas: \$(get_rnas.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}