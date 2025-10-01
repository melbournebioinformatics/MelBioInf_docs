class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - cutadapt

arguments:
  - prefix: '-o'
    valueFrom: out.1.fastq
  - prefix: '-p'
    valueFrom: out.2.fastq

inputs:
  quality_cutoff:
    type: int
    inputBinding:
      prefix: -q

  input_fastq:
    type: File[]
    inputBinding:
      position: 0

outputs:
  trimmed_reads:
    type: File[]
    outputBinding:
      glob: 'out.*.fastq'

requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/cutadapt:1.16--py36_1
