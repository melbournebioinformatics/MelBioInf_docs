class: CommandLineTool

cwlVersion: v1.0

baseCommand:
  - samtools
  - sort

inputs:
  - id: alignment
    type: File
    inputBinding:
      position: 0

outputs:
  - id: sorted_alignment
    type: stdout

requirements:
  - class: DockerRequirement
    dockerPull: biocontainers/samtools
