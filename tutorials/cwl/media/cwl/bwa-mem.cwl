class: CommandLineTool
cwlVersion: v1.0
baseCommand:
  - bwa
  - mem
inputs:
  - id: reads
    type: 'File[]'
    inputBinding:
      position: 1
  - id: reference
    type: File
    inputBinding:
      position: 0
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
outputs:
  - id: alignment
    type: File
    outputBinding:
      glob: output.bam
requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/bwa:0.7.8--hed695b0_5
stdout: output.bam
