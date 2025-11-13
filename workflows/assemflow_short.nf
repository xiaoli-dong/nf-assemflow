/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap } from 'plugin/nf-schema'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

//
// MODULE: nf-core modules
//
include { MULTIQC } from '../modules/local/multiqc/main'

//
//local subworkflows
//
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_assemflow_pipeline'
include { ASSEMBLE_ILLUMINA } from '../subworkflows/local/assembly_illumina'
include { PUBLISH_ASSEMBLIES } from '../modules/local/publish/assemblies'
include { PUBLISH_SAMPLESHEET } from '../modules/local/publish/samplesheet'
include { DEPTH_ILLUMINA } from '../subworkflows/local/depth_illumina'

//
// MODULE: local modules
//
include { CHECKM2_PREDICT } from '../modules/local/checkm2/predict.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMFLOW_SHORT {
    take:
    short_reads

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    short_reads = short_reads.filter { _meta, reads ->
        def readsList = reads instanceof List ? reads : [reads]
        readsList.size() > 0 && readsList.every { it != null && it.exists() && it.size() > 0 }
    }

    // assembly
    ASSEMBLE_ILLUMINA(short_reads)
    // zero size contig can cause some of the program such as bakta, mobsuite file
    ASSEMBLE_ILLUMINA.out.contigs.filter { _meta, contigs -> contigs.countFasta() > 0 }.set { contigs }
    ch_versions = ch_versions.mix(ASSEMBLE_ILLUMINA.out.versions)

    PUBLISH_ASSEMBLIES(contigs)
    assemblies_collected = PUBLISH_ASSEMBLIES.out.contigs.collect().ifEmpty([]).map { it.collate(2) }
    assemblies_collected.view()
    // Generate samplesheet
    PUBLISH_SAMPLESHEET(assemblies_collected)

    //
    // depth report
    //
    if (!params.skip_depth_and_coverage_illumina) {
        short_reads.join(contigs)
            .multiMap { it ->
                reads: [it[0], it[1]]
                contigs: [it[0], it[2]]
            }
            .set {
                ch_input_depth
            }
        DEPTH_ILLUMINA(ch_input_depth.reads, ch_input_depth.contigs)
    }


    if (!params.skip_checkm2) {
        ch_input_checkm2 = contigs.map { _meta, mycontigs -> mycontigs }.collect()
        .map {
            files -> tuple([id: "checkm2"], files)
        }
        //.view()
        CHECKM2_PREDICT(ch_input_checkm2, params.checkm2_db)
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'assemflow_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? channel.fromPath(params.multiqc_config, checkIfExists: true)
        : channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
