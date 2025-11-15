include { FLYE } from '../../modules/local/flye/main'

include { CIRCULARRECENTER_MIDSTARTFLYE } from '../../modules/local/circularrecenter/midstartflye/main'
include { CIRCULARRECENTER_DNAAPLER } from '../../modules/local/circularrecenter/dnaapler/main.nf'
include { GFATOOLS_GFA2FA  } from '../../modules/local/gfatools/gfa2fa/main.nf'
include { POLISHER_NANOPORE } from '../../subworkflows/local/polisher_nanopore'
include { POLISHER_ILLUMINA } from '../../subworkflows/local/polisher_illumina'

include { LRGE } from '../../modules/local/lrge/main'
include { RASUSA_READS } from '../../modules/local/rasusa/reads/main'



workflow ASSEMBLE_HYBRID {
    take:
    long_reads
    short_reads // can be empty for only nanopore data assemflow_long

    main:
    ch_software_versions = channel.empty()
    contigs = channel.empty()
    ASSEMBLE_NANOPORE(long_reads)
    contigs = ASSEMBLE_NANOPORE.out.contigs.filter { meta, mycontigs -> mycontigs.countFasta() > 0 }

    ch_software_versions = ch_software_versions.mix(ASSEMBLE_NANOPORE.out.versions)

    if (params.platform == 'hybrid' & !params.skip_illumina_reads_polish) {
        POLISHER_ILLUMINA(contigs, short_reads)
        ch_software_versions = ch_software_versions.mix(POLISHER_ILLUMINA.out.versions)
        contigs = POLISHER_ILLUMINA.out.contigs

    }

    emit:
    contigs
    versions = ch_software_versions
}



workflow ASSEMBLE_NANOPORE {
    take:
    long_reads

    main:
    ch_software_versions = channel.empty()
    contigs = channel.empty()
    assembly_info = channel.empty()

    //estimate genome size first
    LRGE(long_reads)
    ch_software_versions = ch_software_versions.mix(LRGE.out.versions.first())
    //subsample reads to desired coverage using rasusa to user define coverage (eg. 100x)
    RASUSA_READS(long_reads.join(LRGE.out.gsize))
    long_reads = RASUSA_READS.out.reads
    ch_software_versions = ch_software_versions.mix(RASUSA_READS.out.versions.first())


    //Flye to be the best-performing bacterial genome assembler in many metrics
    if (params.long_read_assembly == 'flye') {

        long_reads.join(LRGE.out.gsize).view()
        FLYE(long_reads.join(LRGE.out.gsize))
        FLYE.out.fasta.filter { _meta, fasta -> fasta.countFasta() > 0 }.set { contigs }
        FLYE.out.txt.filter { _meta, txt -> txt.countLines() > 0 }.set { assembly_info }
        gfa = FLYE.out.gfa

        ch_software_versions = ch_software_versions.mix(FLYE.out.versions.first())

        //recenter the genome before medaka hopes to fix the termial errors
        if (!params.skip_recenter_genome) {
            if (params.recenter_method == 'midstartflye') {
                log.info("About to run midstartflye recentering")
                CIRCULARRECENTER_MIDSTARTFLYE(contigs.join(assembly_info))
                ch_software_versions = ch_software_versions.mix(CIRCULARRECENTER_MIDSTARTFLYE.out.versions.first())
                CIRCULARRECENTER_MIDSTARTFLYE.out.fasta.filter { _meta, fasta -> fasta.countFasta() > 0 }.set { contigs }

            }
            else if (params.recenter_method == 'dnaapler') {
                log.info("About to run dnaapler recentering")
                CIRCULARRECENTER_DNAAPLER(gfa)
                GFATOOLS_GFA2FA(CIRCULARRECENTER_DNAAPLER.out.gfa)
                ch_software_versions = ch_software_versions.mix(CIRCULARRECENTER_DNAAPLER.out.versions.first())
                ch_software_versions = ch_software_versions.mix(GFATOOLS_GFA2FA.out.versions.first())
                GFATOOLS_GFA2FA.out.fasta.filter { _meta, fasta -> fasta.countFasta() > 0 }.set { contigs}

            }
            else {
                log.warn("Unknown recenter_method: ${params.recenter_method}. Skipping recentering step.")
            }
        }
    }

    if (!params.skip_nanopore_reads_polish) {
        log.info("About to run nanopore_polish")
        POLISHER_NANOPORE(contigs, long_reads)
        ch_software_versions = ch_software_versions.mix(POLISHER_NANOPORE.out.versions)
        contigs = POLISHER_NANOPORE.out.contigs

    }

    emit:
    contigs
    versions = ch_software_versions
}
