process RIBOTISH_PREDICT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/../environment.yml"
    container "local.io/ribotish:0.2.7"

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)
    tuple val(meta), path(para)
    path gtf
    path fasta

    output:
    tuple val(meta), path("*.ribotish.txt")     , emit: results
    tuple val(meta), path("*.ribotish_all.txt") , emit: allresults
    tuple val(meta), path(bam)                  , emit: bam
    tuple val(meta), path(bai)                  , emit: bai

    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ribotish predict \\
        -b ${bam} \\
        --ribopara ${para} \\
        -g ${gtf} \\
        -f ${fasta} \\
        -o ${prefix}.ribotish.txt \\
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