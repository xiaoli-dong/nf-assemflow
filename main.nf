#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    xiaoli-dong/nf-assemflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/xiaoli-dong/nf-assemflow
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~_~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ASSEMFLOW_LONG } from './workflows/assemflow_long'
include { ASSEMFLOW_SHORT } from './workflows/assemflow_short'
include { ASSEMFLOW_HYBRID } from './workflows/assemflow_hybrid'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_assemflow_pipeline'
include { PIPELINE_COMPLETION } from './subworkflows/local/utils_nfcore_assemflow_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ABPROVLAB_ASSEMFLOW {
    take:
    short_reads // channel: samplesheet read in from --input
    long_reads // channel:

    main:

    multiqc_report = channel.empty()
    //
    // WORKFLOW: Run pipeline
    //
    if (params.platform == 'illumina') {
        ASSEMFLOW_SHORT(short_reads)
        multiqc_report = ASSEMFLOW_SHORT.out.multiqc_report
    }
    else if (params.platform == 'nanopore') {
        ASSEMFLOW_LONG(long_reads)
        multiqc_report =ASSEMFLOW_LONG.out.multiqc_report
    }
    else if (params.platform == 'hybrid') {
        ASSEMFLOW_HYBRID(short_reads, long_reads)
        multiqc_report = ASSEMFLOW_HYBRID.out.multiqc_report
    }

    emit:
    multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
    )

    //
    // WORKFLOW: Run main workflow
    //
    ABPROVLAB_ASSEMFLOW(
        PIPELINE_INITIALISATION.out.short_reads,
        PIPELINE_INITIALISATION.out.long_reads,
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        ABPROVLAB_ASSEMFLOW.out.multiqc_report,
    )
}
