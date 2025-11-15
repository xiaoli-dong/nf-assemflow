process FLYE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye%3A2.9.6--py39h475c85d_0' :
        'biocontainers/flye%3A2.9.6--py39h475c85d_0' }"

    input:
    tuple val(meta), path(reads), path(gsize)

    output:
    tuple val(meta), path("*.fasta.gz"), emit: fasta
    tuple val(meta), path("*.gfa.gz"), emit: gfa
    tuple val(meta), path("*.gv.gz"), emit: gv
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Run Flye with the converted genome size
    """
    # Read the genome size from the gsize file (base pairs)
    genome_size=\$(cat ${gsize})

    # Convert genome size from base pairs (bp) to megabases (Mb)
    genome_size_in_m=\$((genome_size / 1000000))  # Convert to Mb

    # Use megabases (Mb) for Flye (uncomment if you want to use Gb instead)
    final_genome_size=\${genome_size_in_m}m

    flye $args \\
        $reads \\
        --out-dir . \\
        --threads $task.cpus \\
        --genome-size \${final_genome_size} \\
        $args2 > flye.log 2>&1

    gzip -c assembly.fasta > ${prefix}.assembly.fasta.gz
    gzip -c assembly_graph.gfa > ${prefix}.assembly_graph.gfa.gz
    gzip -c assembly_graph.gv > ${prefix}.assembly_graph.gv.gz
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv flye.log ${prefix}.flye.log
    mv params.json ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub | gzip -c > ${prefix}.assembly.fasta.gz
    echo stub | gzip -c > ${prefix}.assembly_graph.gfa.gz
    echo stub | gzip -c > ${prefix}.assembly_graph.gv.gz
    echo contig_1 > ${prefix}.assembly_info.txt
    echo stub > ${prefix}.flye.log
    echo stub > ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version)
    END_VERSIONS
    """
}
