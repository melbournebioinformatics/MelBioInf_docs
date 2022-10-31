class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
baseCommand:
  - freebayes
inputs:
  - id: reference
    type: File
    inputBinding:
      position: 0
      prefix: '--fasta-reference'
  - id: bam
    type: File
    inputBinding:
      position: 0
outputs:
  - id: output
    type: File
    outputBinding:
      glob: output.vcf
requirements:
  - class: DockerRequirement
    dockerPull: maxulysse/freebayes
stdout: output.vcf
