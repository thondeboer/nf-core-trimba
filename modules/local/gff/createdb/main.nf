process GFF_CREATEDB {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffutils:0.11.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/gffutils:0.12--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.db"), emit: db

    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gff_createdb.py \\
        --gtf $gtf \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcreatedb: \$(gff_createdb.py --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcreatedb: \$(gff_createdb.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}