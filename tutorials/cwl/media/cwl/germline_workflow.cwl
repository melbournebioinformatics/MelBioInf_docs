class: Workflow
cwlVersion: v1.0

inputs:
  - id: reference
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa

  - id: reads
    type: 'File[]'

outputs:
  - id: output
    outputSource:
      - freebayes/output
    type: File

steps:
  - id: bwa_mem
    in:
      - id: reads
        source:
          - reads
      - id: reference
        source: reference
    out:
      - id: alignment
    run: ./bwa-mem.cwl

  - id: samtools_sort
    in:
      - id: alignment
        source: bwa_mem/alignment
    out:
      - id: sorted_alignment
    run: ./samtools-sort.cwl

  - id: samtools_index
    in:
      - id: alignment
        source: samtools_sort/sorted_alignment
    out:
      - id: alignment_with_index
    run: ./samtools-index.cwl

  - id: freebayes
    in:
      - id: reference
        source: reference
      - id: bam
        source: samtools_index/alignment_with_index
    out:
      - id: output
    run: ./freebayes.cwl
