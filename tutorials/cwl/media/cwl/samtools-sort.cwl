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
    type: File
    outputBinding:
      glob: sorted_alignment.bam

stdout: sorted_alignment.bam

requirements:
  - class: DockerRequirement
    dockerPull: biocontainers/samtools
