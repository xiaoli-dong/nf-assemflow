process BBMAP_BBNORM {
    tag "$meta.id"
    label 'process_medium'

     conda "${moduleDir}/environment.yml"
     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bbmap%3A39.37--he5f24ec_0':
        'biocontainers/bbmap%3A39.37--he5f24ec_0' }"


    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    input  = meta.single_end ? "in=${fastq.join(',')}" : "in=${fastq[0]} in2=${fastq[1]}"
    output = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_1.nm.fastq.gz out2=${prefix}_2.nm.fastq.gz"

    memory = '-Xmx3g'
    if ( ! task.memory ) {
        log.info '[BBNorm]: Available memory not known, defaulting to 3 GB. Specify process memory requirements to change this.'
    } else {
        memory = "-Xmx${Math.round(Math.max(1, Math.floor(task.memory.toGiga() * 0.95)))}g"
    }

    """
    bbnorm.sh \\
        $input \\
        $output \\
        hist=${prefix}.hist \\
        $args \\
        threads=$task.cpus \\
        $memory \\
        &> ${prefix}.bbnorm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
