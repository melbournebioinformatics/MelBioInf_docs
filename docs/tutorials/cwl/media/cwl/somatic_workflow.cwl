class: Workflow
cwlVersion: v1.0

inputs:
  - id: normal_reads
    type: 'File[]'

  - id: reference
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  - id: tumour_reads
    type: 'File[]'

outputs:
  - id: vcf
    outputSource:
      - somatic_sniper/vcf
    type: File

steps:
  - id: normal_alignment
    in:
      - id: reads
        source:
          - normal_reads
      - id: reference
        source: reference
    out:
      - id: alignment_with_index
    run: alignment_workflow.cwl

  - id: tumour_alignment
    in:
      - id: reads
        source:
          - tumour_reads
      - id: reference
        source: reference
    out:
      - id: alignment_with_index
    run: alignment_workflow.cwl

  - id: somatic_sniper
    in:
      - id: normal
        source: normal_alignment/alignment_with_index
      - id: reference
        source: reference
      - id: tumour
        source: tumour_alignment/alignment_with_index
    out:
      - id: vcf
    run: ./somatic-sniper.cwl

requirements:
  - class: SubworkflowFeatureRequirement
