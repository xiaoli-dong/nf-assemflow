include { ASSEMBLYSTATS } from '../../modules/local/stats/assemblystats'
include { REFORMATASSEMBLYSTATS } from '../../modules/local/stats/reformatassemblystats'

workflow POLISHER_ILLUMINA {
    take:
    draft_contigs
    short_reads

    main:

    ch_software_versions = channel.empty()
    contigs = draft_contigs

    if (!params.skip_polypolish) {
        RUN_POLYPOLISH(draft_contigs, short_reads)
        contigs = draft_contigs.merge(RUN_POLYPOLISH.out.contigs.ifEmpty([]))
            .map { it ->
                if (it.size() == 4) {
                    [it[2], it[3]]
                }
                else {
                    [it[0], it[1]]
                }
            }

        ch_software_versions = ch_software_versions.mix(RUN_POLYPOLISH.out.versions)
    }

    if (!params.skip_pypolca) {

        RUN_PYPOLCA(contigs, short_reads)
        contigs = contigs
            .merge(RUN_PYPOLCA.out.contigs.ifEmpty([]))
            .map { it ->
                if (it.size() == 4) {
                    [it[2], it[3]]
                }
                else {
                    [it[0], it[1]]
                }
            }

        ch_software_versions = ch_software_versions.mix(RUN_PYPOLCA.out.versions)
    }

    ASSEMBLYSTATS(contigs)
    REFORMATASSEMBLYSTATS(ASSEMBLYSTATS.out.stats)
    stats = REFORMATASSEMBLYSTATS.out.tsv

    emit:
    contigs
    stats
    versions = ch_software_versions
}


include { BWAMEM2_INDEX } from '../../modules/local/bwamem2/index/main'
include {
    BWAMEM2_MEM as BWAMEM2_MEM_1 ;
    BWAMEM2_MEM as BWAMEM2_MEM_2
} from '../../modules/local/bwamem2/mem'
include { POLYPOLISH } from '../../modules/local/polypolish'

workflow RUN_POLYPOLISH {
    take:
    draft_contigs
    short_reads

    main:
    ch_software_versions = channel.empty()

    //draft_contigs.view()
    short_reads
        .map { meta, reads ->
            def new_meta = [:]
            new_meta.id = meta.id
            [new_meta, meta, reads]
        }
        .set {
            ch_input_reads
        }

    draft_contigs
        .map { meta, contigs ->
            def new_meta = [:]
            new_meta.id = meta.id
            [new_meta, meta, contigs]
        }
        .set {
            ch_input_contigs
        }

    ch_input_reads
        .join(ch_input_contigs)
        .multiMap { it ->
            reads: [it[1], it[2]]
            draft_contigs: [it[1], it[4]]
        }
        .set {
            ch_input_all
        }

    BWAMEM2_INDEX(ch_input_all.draft_contigs)
    ch_software_versions = ch_software_versions.mix(BWAMEM2_INDEX.out.versions.first())

    ch_input_all.reads
        .join(BWAMEM2_INDEX.out.index)
        .multiMap { it ->
            read_1: [it[0], it[1][0]]
            read_2: [it[0], it[1][1]]
            bwa_index: [it[0], it[2]]
        }
        .set {
            ch_input
        }

    BWAMEM2_MEM_1(ch_input.read_1, ch_input.bwa_index, false)
    BWAMEM2_MEM_2(ch_input.read_2, ch_input.bwa_index, false)
    ch_software_versions = ch_software_versions.mix(BWAMEM2_MEM_1.out.versions.first())

    ch_input_all.draft_contigs
        .join(BWAMEM2_MEM_1.out.sam)
        .join(BWAMEM2_MEM_2.out.sam)
        .multiMap { it ->
            contigs: [it[0], it[1]]
            sam: [it[0], it[2], it[3]]
        }
        .set {
            ch_input
        }
    POLYPOLISH(ch_input.contigs, ch_input.sam)
    ch_software_versions = ch_software_versions.mix(POLYPOLISH.out.versions.first())

    contigs = POLYPOLISH.out.contigs
    ASSEMBLYSTATS(contigs)
    REFORMATASSEMBLYSTATS(ASSEMBLYSTATS.out.stats)
    stats = REFORMATASSEMBLYSTATS.out.tsv

    emit:
    contigs
    stats
    versions = ch_software_versions
}

include { PYPOLCA } from '../../modules/local/pypolca'

workflow RUN_PYPOLCA {
    take:
    draft_contigs
    short_reads

    main:
    ch_software_versions = channel.empty()

    //draft_contigs.view()
    short_reads
        .map { meta, reads ->
            def new_meta = [:]
            new_meta.id = meta.id
            [new_meta, meta, reads]
        }
        .set {
            ch_input_reads
        }

    draft_contigs
        .map { meta, contigs ->
            def new_meta = [:]
            new_meta.id = meta.id
            [new_meta, meta, contigs]
        }
        .set {
            ch_input_contigs
        }

    ch_input_reads
        .join(ch_input_contigs)
        .multiMap { it ->
            reads: [it[1], it[2]]
            draft_contigs: [it[1], it[4]]
        }
        .set {
            ch_input
        }


    PYPOLCA(ch_input.reads, ch_input.draft_contigs)
    ch_software_versions = ch_software_versions.mix(PYPOLCA.out.versions.first())

    contigs = PYPOLCA.out.corrected_contigs
    ASSEMBLYSTATS(contigs)
    REFORMATASSEMBLYSTATS(ASSEMBLYSTATS.out.stats)
    stats = REFORMATASSEMBLYSTATS.out.tsv


    emit:
    contigs
    stats
    versions = ch_software_versions
}
