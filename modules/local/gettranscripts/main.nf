process GETTRANSCRIPTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/gffutils"

    input:
    tuple val(meta), path(ribotish)
    tuple val(meta), path(bedgraph)
    path(dbfile)
    path(fasta)
    path(genes)

    output:
    tuple val(meta), path("*.fa")                  , emit: fasta   , optional: true
    tuple val(meta), path("*.bedgraph")            , emit: bedgraph, optional: true
    tuple val(meta), path("*filtered.ribotish.txt"), emit: ribotish, optional: true

    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_transcripts.py \\
        --db_file ${dbfile} \\
        --fasta ${fasta} \\
        --genes_file ${genes} \\
        --ribotish_file ${ribotish} \\
        --bedgraph ${bedgraph} \\
        --prefix ${prefix}.fiveprimeutrs \\
        ${args}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gettranscripts: \$(get_transcripts.py --version | cut -d' ' -f2)
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gettranscripts: \$(get_transcripts.py --version | cut -d' ' -f2)
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

