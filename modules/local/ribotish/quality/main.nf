process RIBOTISH_QUALITY {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ribotish:0.2.5--pyh864c0ab_1' :
        'quay.io/biocontainers/ribotish:0.2.7--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    path(gtf)

    output:
    tuple val(meta), path("*.txt"), emit: results, optional: true
    tuple val(meta), path("*.pdf"), emit: pdf, optional: true
    tuple val(meta), path("*.py") , emit: psites, optional: true

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ribotish quality \\
        -b ${bam} \\
        -g ${gtf} \\
        -o ${prefix}_qual.txt \\
        -f ${prefix}_qual.pdf \\
        -r ${prefix}_para.py \\
        -p $task.cpus \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ribotish: \$(ribotish --version | cut -d' ' -f2)
    END_VERSIONS
    """
}