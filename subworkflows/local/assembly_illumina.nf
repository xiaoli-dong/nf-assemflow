//https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8380430/
//https://peerj.com/articles/12446.pdf

include { SPADES } from '../../modules/local/spades/main'
include { UNICYCLER } from '../../modules/local/unicycler/main'
include { MEGAHIT } from '../../modules/local/megahit/main'
include { SKESA } from '../../modules/local/skesa'
include { SHOVILL } from '../../modules/local/shovill/main'
include { BBMAP_BBNORM } from '../../modules/local/bbmap/bbnorm/main'

workflow ASSEMBLE_ILLUMINA {
    take:
    reads

    main:
    ch_versions = channel.empty()
    contigs = channel.empty()


    if (params.short_read_assembly != 'shovill') {
        //https://help.geneious.com/hc/en-us/articles/360044626852-Best-practice-for-preprocessing-NGS-reads-in-Geneious-Prime
        BBMAP_BBNORM(reads)
        reads = BBMAP_BBNORM.out.fastq

        if (params.short_read_assembly == 'spades') {
            reads
                .map { meta, myreads -> [meta, myreads, [], []] }
                .set { input }
            SPADES(input, [], [])

            //get rid of zero size contig file and avoid the downstream crash
            SPADES.out.contigs
                .filter { meta, mycontigs -> mycontigs.countFasta() > 0 }
                .set { contigs }

            ch_versions = ch_versions.mix(SPADES.out.versions.first())
        }
        else if (params.short_read_assembly == 'skesa') {
            SKESA(reads)
            SKESA.out.contigs.filter { meta, mycontigs -> mycontigs.countFasta() > 0 }.set { contigs }
            ch_versions = ch_versions.mix(SKESA.out.versions.first())
        }
        else if (params.short_read_assembly == 'unicycler') {
            reads.map { meta, myreads -> [meta, myreads, []] }.set { input }
            UNICYCLER(input)
            UNICYCLER.out.scaffolds.filter { meta, scaffolds -> scaffolds.size() > 0 }.set { contigs }
            ch_versions = ch_versions.mix(UNICYCLER.out.versions.first())
        }
        else if (params.short_read_assembly == 'megahit') {
            MEGAHIT(reads)
            MEGAHIT.out.contigs.filter { meta, mycontigs -> mycontigs.countFasta() > 0 }.set { contigs }
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        }
    }
    else if (params.short_read_assembly == 'shovill') {
        SHOVILL(reads)
        SHOVILL.out.contigs.filter { meta, mycontigs -> mycontigs.countFasta() > 0 }.set { contigs }
        ch_versions = ch_versions.mix(SHOVILL.out.versions.first())
    }

    emit:
    contigs
    versions = ch_versions
}
