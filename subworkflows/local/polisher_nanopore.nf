include { ASSEMBLYSTATS } from '../../modules/local/stats/assemblystats'
include { REFORMATASSEMBLYSTATS } from '../../modules/local/stats/reformatassemblystats'
include { MEDAKA_CONSENSUS } from '../../modules/local/medaka/consensus'
include { DORADO_ALIGNER } from '../../modules/local/dorado/aligner'
include { DORADO_POLISH } from '../../modules/local/dorado/polish'

workflow POLISHER_NANOPORE {
    take:
    draft_contigs
    long_reads

    main:
    ch_software_versions = channel.empty()
    contigs = draft_contigs

    if(params.nanopore_reads_polisher == 'medaka'){
        MEDAKA_CONSENSUS(long_reads.join(draft_contigs))
        ch_software_versions = ch_software_versions.mix(MEDAKA_CONSENSUS.out.versions.first())
        contigs = MEDAKA_CONSENSUS.out.assembly
    }
    //todo
    /*
    too buggy yet and not implement it for now
    https://github.com/nanoporetech/dorado/issues/1384
    */
    /* else if(params.nanopore_reads_polisher == 'dorado'){
        DORADO_ALIGNER(long_reads.join(draft_contigs))
        ch_software_versions = ch_software_versions.mix(DORADO_ALIGNER.out.versions.first())

        DORADO_POLISH(DORADO_ALIGNER.out.bam.join(DORADO_ALIGNER.out.bai).join(draft_contigs))
        ch_software_versions = ch_software_versions.mix(DORADO_POLISH.out.versions.first())
        contigs = DORADO_POLISH.out.assembly
    }
 */
    ASSEMBLYSTATS(contigs)
    ch_software_versions = ch_software_versions.mix(ASSEMBLYSTATS.out.versions.first())

    REFORMATASSEMBLYSTATS (ASSEMBLYSTATS.out.stats)
    ch_software_versions = ch_software_versions.mix(REFORMATASSEMBLYSTATS.out.versions.first())
    stats = REFORMATASSEMBLYSTATS.out.tsv

    emit:
    contigs
    stats
    versions = ch_software_versions

}

