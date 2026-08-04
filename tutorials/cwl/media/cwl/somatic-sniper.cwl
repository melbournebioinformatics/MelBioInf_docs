class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - bam-somaticsniper

inputs:
  - id: normal
    type: File
    inputBinding:
      position: 1
    label: Normal BAM

  - id: reference
    type: File
    inputBinding:
      position: 0
      prefix: '-f'
    label: Reference genome

  - id: tumour
    type: File
    inputBinding:
      position: 0
    label: Tumour BAM

outputs:
  - id: vcf
    type: stdout

arguments:
  - position: 2
    valueFrom: output.vcf

requirements:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/somatic-sniper:1.0.5.0--0'

stdout: variants.vcf
