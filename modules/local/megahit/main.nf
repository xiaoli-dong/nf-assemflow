process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/megahit%3A1.2.9--haf24da9_8'
        : 'biocontainers/megahit%3A1.2.9--haf24da9_8'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa.gz"), emit: contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa.gz"), emit: k_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa.gz"), emit: addi_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa.gz"), emit: local_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        megahit \\
            -r ${reads} \\
            -t ${task.cpus} \\
            ${args} \\
            -o megahit_out \\
            --out-prefix ${prefix}

        gzip --no-name ${args2} megahit_out/*.fa megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
    else {
        """
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -t ${task.cpus} \\
             ${args} \\
            -o megahit_out \\
            --out-prefix ${prefix}

        gzip --no-name ${args2} megahit_out/*.fa megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }
}
