class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - cutadapt

inputs:
  quality_cutoff:
    type: int
    inputBinding:
      prefix: -q

  input_fastq:
    type: File
    inputBinding:
      position: 0

outputs:
  trimmed_reads:
    type: stdout

requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/cutadapt:1.16--py36_1
