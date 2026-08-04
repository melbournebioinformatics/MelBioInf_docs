class: Workflow

cwlVersion: v1.0

inputs:
  reference:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  reads:
    type: 'File[]'

outputs:
  alignment_with_index:
    outputSource:
      - samtools_index/alignment_with_index
    type: File

steps:
  - id: bwa_mem
    in:
      reads: reads
      reference: reference
    out:
      - alignment
    run: bwa-mem.cwl

  - id: samtools_sort
    in:
      alignment: bwa_mem/alignment
    out:
      - sorted_alignment
    run: samtools-sort.cwl

  - id: samtools_index
    in:
      alignment: samtools_sort/sorted_alignment
    out:
      - alignment_with_index
    run: samtools-index.cwl
