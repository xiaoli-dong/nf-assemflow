process RASUSA_READS {
    tag "$meta.id"
    label 'process_medium'


    conda "bioconda::rasusa=2.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rasusa:2.2.2--hc1c3326_0':
        'biocontainers/rasusa:2.2.2--hc1c3326_0' }"

    input:
    tuple val(meta), path(reads), path(gsize)


    output:
    tuple val(meta), path("*rasusa_subsample.fastq.gz"), emit: reads
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    //subsample
    // cannot use gzip contig file as  for polca, otherwise it will hang and never finish
    """
    rasusa reads \\
        $args \\
        --genome-size \$(cat ${gsize}) \\
        -O g \\
        -o ${prefix}.rasusa_subsample.fastq.gz \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(echo \$(rasusa -V 2>&1) | sed 's/^rasusa //;')
    END_VERSIONS
    """
}
