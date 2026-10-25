//
// Subworkflow with functionality specific to the nf-core/pdefinisher pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Display help and exit if requested
    //
    if (help || help_full) {
        log.info paramsHelp(monochrome_logs)
        exit 0
    }

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Pre-run verification gate: fail-fast parameter & path integrity check
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, assembly, bam, bai ->
            def resolved_bai = bai ? file(bai, checkIfExists: true) : (
                file("${bam}.bai").exists() ? file("${bam}.bai", checkIfExists: true) :
                file("${bam}.csi", checkIfExists: true)
            )
            def norm_chroms = meta.chromosomes ? meta.chromosomes.replaceAll(';', ',') : (params.chromosomes ?: null)
            def norm_telo   = meta.telomere_repeat ?: (params.telomere_repeat ?: 'TTAGGC')
            def enriched_meta = meta + [
                chromosomes    : norm_chroms,
                telomere_repeat: norm_telo
            ]
            return [ enriched_meta, file(assembly, checkIfExists: true), file(bam, checkIfExists: true), resolved_bai ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow)
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Command-line help banner
//
def paramsHelp(monochrome_logs=false) {
    def logo = """
------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\\ | |__  __ /  ` /  \\ |__) |__         }  {
  | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                        `._,._,'
  nf-core/pdefinisher ${workflow.manifest.version ?: '1.0.0'}
------------------------------------------------------
"""
    def help_text = """
${logo}
Automated T2T and PDE finishing for genomes undergoing Programmed DNA Elimination

Usage:
  nextflow run nf-core/pdefinisher --input <samplesheet.csv> --outdir <OUTDIR> [options]

Mandatory arguments:
  --input                Path to comma-separated samplesheet CSV
  --outdir               Path to directory where results are published

PDE Finisher options:
  --telomere_repeat      Canonical telomere repeat monomer (default: 'TTAGGC')
  --chromosomes          Comma-separated list of chromosome names (e.g. 'I,II,III,IV,V,X')
  --curation_policy      Path to curation policy JSON (default: 'assets/curation_policy.json')
  --skip_validation      Skip BAM re-mapping and Sniffles SV validation (default: true)
  --run_motif_discovery  Run MEME grid search and FIMO whole-genome scanning (default: false)
  --fimo_pvalue_threshold P-value threshold for FIMO scanning (default: 1e-4)

Tool-specific options:
  --delfies_args         Extra arguments for DELFIES breakpoint detection
  --mosdepth_args        Extra arguments for mosdepth coverage profiling
  --hifiasm_args         Extra arguments for hifiasm assembly
  --minimap2_align_args  Extra arguments for minimap2 alignment

Generic options:
  --help                 Display this help message and exit
  --monochrome_logs      Disable ANSI colour formatting in logs
===================================================================
""".stripIndent()

    return help_text
}

//
// Native samplesheet CSV parser and validator (runs without nf-schema plugin)
//
def samplesheetToList(String samplesheet_file, String schema_path = null) {
    def ss = file(samplesheet_file)
    def lines = ss.readLines().findAll { it.trim() && !it.trim().startsWith('#') }
    def separator = samplesheet_file.endsWith('.tsv') ? '\t' : ','
    def header = lines[0].tokenize(separator).collect { it.trim() }
    def sample_idx      = header.indexOf('sample')
    def assembly_idx    = header.indexOf('assembly')
    def bam_idx         = header.indexOf('bam')
    def bai_idx         = header.indexOf('bai')
    def chromosomes_idx = header.indexOf('chromosomes')
    def telo_idx        = header.indexOf('telomere_repeat')

    def records = []
    lines.drop(1).each { line ->
        def cols = line.tokenize(separator).collect { it.trim() }
        def sample = sample_idx >= 0 && sample_idx < cols.size() ? cols[sample_idx] : null
        def assembly = assembly_idx >= 0 && assembly_idx < cols.size() ? cols[assembly_idx] : null
        def bam = bam_idx >= 0 && bam_idx < cols.size() ? cols[bam_idx] : null
        def bai = bai_idx >= 0 && bai_idx < cols.size() && cols[bai_idx] ? cols[bai_idx] : null
        def chroms = chromosomes_idx >= 0 && chromosomes_idx < cols.size() && cols[chromosomes_idx] ? cols[chromosomes_idx] : null
        def telo = telo_idx >= 0 && telo_idx < cols.size() && cols[telo_idx] ? cols[telo_idx] : null

        def meta = [ id: sample ]
        if (chroms) meta.chromosomes = chroms
        if (telo)   meta.telomere_repeat = telo

        records << [ meta, assembly, bam, bai ]
    }
    return records
}

//
// Native parameter summary map for MultiQC and completion reports
//
def paramsSummaryMap(workflow, Map options = [:]) {
    def summary = [:]
    summary['Input / Output Options'] = [
        'input'                      : params.input,
        'outdir'                     : params.outdir,
        'email'                      : params.email,
        'multiqc_title'              : params.multiqc_title
    ]
    summary['PDE Finisher Options'] = [
        'telomere_repeat'            : params.telomere_repeat,
        'chromosomes'                : params.chromosomes,
        'curation_policy'            : params.curation_policy,
        'skip_validation'            : params.skip_validation,
        'run_motif_discovery'        : params.run_motif_discovery,
        'fimo_pvalue_threshold'      : params.fimo_pvalue_threshold
    ]
    summary['Tool Arguments'] = [
        'delfies_args'               : params.delfies_args,
        'mosdepth_args'              : params.mosdepth_args,
        'hifiasm_args'               : params.hifiasm_args,
        'minimap2_align_args'        : params.minimap2_align_args
    ]
    summary['Execution Environment'] = [
        'launchDir'                  : workflow.launchDir,
        'workDir'                    : workflow.workDir,
        'projectDir'                 : workflow.projectDir,
        'profile'                    : workflow.profile,
        'containerEngine'            : workflow.containerEngine
    ]
    return summary
}

//
// Pre-run verification gate: fail-fast parameter & path integrity check
//
def validateInputParameters() {
    def validation_errors = []

    // 1. Mandatory Input Samplesheet
    if (!params.input) {
        validation_errors.add("Mandatory parameter --input <path/to/samplesheet.csv> is missing.")
    } else {
        def ss_file = file(params.input)
        if (!ss_file.exists()) {
            validation_errors.add("Input samplesheet file does not exist: ${params.input}")
        } else if (ss_file.isEmpty()) {
            validation_errors.add("Input samplesheet file is empty: ${params.input}")
        } else {
            def lines = ss_file.readLines().findAll { it.trim() && !it.trim().startsWith('#') }
            if (lines.size() < 2) {
                validation_errors.add("Input samplesheet contains no sample rows: ${params.input}")
            } else {
                def separator = params.input.endsWith('.tsv') ? '\t' : ','
                def header = lines[0].tokenize(separator).collect { it.trim() }
                def sample_idx   = header.indexOf('sample')
                def assembly_idx = header.indexOf('assembly')
                def bam_idx      = header.indexOf('bam')
                def bai_idx      = header.indexOf('bai')
                def telo_idx     = header.indexOf('telomere_repeat')

                if (sample_idx < 0) {
                    validation_errors.add("Samplesheet header is missing required column 'sample': ${lines[0]}")
                }
                if (assembly_idx < 0) {
                    validation_errors.add("Samplesheet header is missing required column 'assembly': ${lines[0]}")
                }
                if (bam_idx < 0) {
                    validation_errors.add("Samplesheet header is missing required column 'bam': ${lines[0]}")
                }

                lines.drop(1).eachWithIndex { line, line_idx ->
                    def row_num = line_idx + 2
                    def cols = line.tokenize(separator).collect { it.trim() }
                    def sample_id = sample_idx >= 0 && sample_idx < cols.size() ? cols[sample_idx] : "row_${row_num}"

                    if (!sample_id || sample_id.contains(' ')) {
                        validation_errors.add("Samplesheet row ${row_num}: 'sample' ID must not be empty or contain whitespace (found: '${sample_id}').")
                    }

                    if (assembly_idx >= 0 && assembly_idx < cols.size()) {
                        def fa_path = cols[assembly_idx]
                        if (!fa_path) {
                            validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): 'assembly' column is empty.")
                        } else {
                            def f = file(fa_path)
                            if (!f.exists()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): assembly file does not exist: ${fa_path}")
                            } else if (f.isEmpty()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): assembly file is empty: ${fa_path}")
                            }
                        }
                    }

                    if (bam_idx >= 0 && bam_idx < cols.size()) {
                        def bam_path = cols[bam_idx]
                        if (!bam_path) {
                            validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): 'bam' column is empty.")
                        } else {
                            def f = file(bam_path)
                            if (!f.exists()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): BAM file does not exist: ${bam_path}")
                            } else if (f.isEmpty()) {
                                validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): BAM file is empty: ${bam_path}")
                            } else {
                                def bai_path = bai_idx >= 0 && bai_idx < cols.size() && cols[bai_idx] ? cols[bai_idx] : null
                                if (bai_path) {
                                    def b = file(bai_path)
                                    if (!b.exists()) {
                                        validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): BAI file does not exist: ${bai_path}")
                                    }
                                } else {
                                    def default_bai = file("${bam_path}.bai")
                                    def default_csi = file("${bam_path}.csi")
                                    if (!default_bai.exists() && !default_csi.exists()) {
                                        validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): Neither ${bam_path}.bai nor ${bam_path}.csi index found.")
                                    }
                                }
                            }
                        }
                    }

                    if (telo_idx >= 0 && telo_idx < cols.size() && cols[telo_idx]) {
                        def telo = cols[telo_idx]
                        if (!telo.matches('^[ACGTacgt]+$')) {
                            validation_errors.add("Samplesheet row ${row_num} (sample '${sample_id}'): telomere_repeat '${telo}' contains invalid characters (must be A, C, G, T).")
                        }
                    }
                }
            }
        }
    }

    // 2. Mandatory Output Directory
    if (!params.outdir) {
        validation_errors.add("Mandatory parameter --outdir <path/to/output_dir> is missing.")
    }

    // 3. Curation Policy verification
    if (params.curation_policy) {
        def cp = file(params.curation_policy)
        if (!cp.exists()) {
            validation_errors.add("Curation policy JSON does not exist: ${params.curation_policy}")
        }
    }

    // Fail-fast error logging
    if (validation_errors) {
        log.error "==================================================================="
        log.error " PRE-RUN VERIFICATION GATE FAILED: ${validation_errors.size()} CRITICAL CONFIGURATION ISSUE(S)"
        log.error "==================================================================="
        validation_errors.eachWithIndex { err, i ->
            log.error " [${i+1}] ${err}"
        }
        log.error "==================================================================="
        error "Pre-run verification gate failed. Please correct the parameters above before executing."
    }

    log.info "==================================================================="
    log.info " PRE-RUN VERIFICATION GATE PASSED: All inputs, samplesheet records, and paths verified."
    log.info "==================================================================="
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
