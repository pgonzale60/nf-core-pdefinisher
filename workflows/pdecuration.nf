/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { AUTOMATED_PDE_CURATION } from '../subworkflows/local/automated_pde_curation/main'
include { CURATIONSUMMARY        } from '../modules/local/curationsummary/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pdecuration_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PDECURATION {

    take:
    ch_samplesheet // channel: [ val(meta), path(assembly), path(bam), path(bai) ]
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    ch_assembly = ch_samplesheet.map { meta, assembly, _bam, _bai -> tuple(meta, assembly) }
    ch_bam      = ch_samplesheet.map { meta, _assembly, bam, bai  -> tuple(meta, bam, bai) }

    //
    // SUBWORKFLOW: Run end-to-end Automated PDE Curation
    //
    AUTOMATED_PDE_CURATION (
        ch_assembly,
        ch_bam,
        params.telomere_repeat ?: 'TTAGGC',
        params.chromosomes ?: ''
    )
    ch_versions = ch_versions.mix(AUTOMATED_PDE_CURATION.out.versions)

    //
    // MODULE: Cohort-level Curation Summary across all processed samples
    //
    ch_all_scorecards  = AUTOMATED_PDE_CURATION.out.scorecard_json.map { _meta, f -> f }.collect()
    ch_all_validations = AUTOMATED_PDE_CURATION.out.validation_json.map { _meta, f -> f }.collect()

    ch_summary_input = ch_all_scorecards
        .combine(ch_all_validations)
        .map { scorecards, validations ->
            tuple([ id: 'cohort' ], scorecards, validations)
        }

    CURATIONSUMMARY (
        ch_summary_input
    )
    ch_versions = ch_versions.mix(CURATIONSUMMARY.out.versions_summary)

    //
    // Collate and save software versions
    //
    def all_raw_versions = ch_versions
        .mix(channel.topic("versions"))
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = all_raw_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(all_raw_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'pdecuration_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'pdecuration'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    curated_fasta   = AUTOMATED_PDE_CURATION.out.curated_fasta
    curated_agp     = AUTOMATED_PDE_CURATION.out.curated_agp
    blueprint_json  = AUTOMATED_PDE_CURATION.out.blueprint_json
    scorecard_json  = AUTOMATED_PDE_CURATION.out.scorecard_json
    validation_json = AUTOMATED_PDE_CURATION.out.validation_json
    summary_md      = CURATIONSUMMARY.out.summary_md
    multiqc_report  = MULTIQC.out.report.map { _meta, report -> [report] }.toList()
    versions        = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
