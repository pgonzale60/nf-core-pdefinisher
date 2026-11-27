include { optionalEvidence } from '../main'

workflow OPTIONAL_EVIDENCE {

    take:
    ch_assembly
    ch_evidence

    main:
    retained = optionalEvidence(ch_assembly, ch_evidence)

    emit:
    samples = retained
}
