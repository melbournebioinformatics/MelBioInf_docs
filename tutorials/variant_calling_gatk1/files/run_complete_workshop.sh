#!/bin/bash

mkdir  "MYTEST_"$RUNNUMBER

cd "MYTEST_"$RUNNUMBER

mkdir data
mkdir output
mkdir reference
mkdir reference/hg38
mkdir scripts
mkdir slurm_scripts
mkdir temp
mkdir tools

cp -p /mnt/shared_data/NA12878.chr20.region_1.fastq.gz data/.
cp -p /mnt/shared_data/NA12878.chr20.region_2.fastq.gz data/.

ln -s /mnt/shared_data/* reference/hg38/.

bwa mem -M -t 2 \
-R "@RG\tID:SRR622461.7\tSM:NA12878\tLB:ERR194147\tPL:ILLUMINA" \
reference/hg38/Homo_sapiens_assembly38.fasta \
data/NA12878.chr20.region_1.fastq.gz \
data/NA12878.chr20.region_2.fastq.gz | \
samtools view -b -h -o output/NA12878.bam -

picard -Xmx7g SortSam \
    I=output/NA12878.bam \
    O=output/NA12878.sort.bam \
    VALIDATION_STRINGENCY=LENIENT \
    SORT_ORDER=coordinate \
    MAX_RECORDS_IN_RAM=3000000 \
    CREATE_INDEX=True

picard -Xmx7g MarkDuplicates \
    I=output/NA12878.sort.bam \
    O=output/NA12878.sort.dup.bam \
    METRICS_FILE=output/marked_dup_metrics.txt

gatk --java-options "-Xmx7g" BaseRecalibrator \
    -I output/NA12878.sort.dup.bam \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    --known-sites reference/hg38/dbsnp_146.hg38.vcf.gz \
    -O output/recal_data.table

gatk --java-options "-Xmx7g" ApplyBQSR \
    -I output/NA12878.sort.dup.bam \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file output/recal_data.table \
    -O output/NA12878.sort.dup.bqsr.bam


fastqc data/NA12878.chr20.region_1.fastq.gz data/NA12878.chr20.region_2.fastq.gz -o output/

picard CollectMultipleMetrics R=reference/hg38/Homo_sapiens_assembly38.fasta I=output/NA12878.sort.dup.bqsr.bam O=output/NA12878.sort.dup.bqsr.CollectMultipleMetrics

multiqc output/. -o output/.


gatk --java-options "-Xmx7g" HaplotypeCaller \
    -I output/NA12878.sort.dup.bqsr.bam \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -ERC GVCF \
    -L chr20 \
    -O output/NA12878.g.vcf.gz

cp /mnt/shared_data/NA12891.g.vcf.gz* output/.

cp /mnt/shared_data/NA12892.g.vcf.gz* output/.

gatk --java-options "-Xmx7g" CombineGVCFs \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/NA12878.g.vcf.gz \
    -V output/NA12891.g.vcf.gz \
    -V output/NA12892.g.vcf.gz \
    -L chr20 \
    -O output/cohort.g.vcf.gz

gatk --java-options "-Xmx7g" GenotypeGVCFs \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/cohort.g.vcf.gz \
    -L chr20 \
    -O output/output.vcf.gz

gatk --java-options "-Xmx7g" VariantRecalibrator \
    -V output/output.vcf.gz \
    --trust-all-polymorphic \
    -mode SNP \
    --max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 reference/hg38/hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=true,prior=12 reference/hg38/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10 reference/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7 reference/hg38/dbsnp_138.hg38.vcf.gz \
    -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
    -O output/cohort_snps.recal \
    --tranches-file output/cohort_snps.tranches

gatk --java-options "-Xmx7g" ApplyVQSR \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/output.vcf.gz \
    -O output/output.vqsr.vcf \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file output/cohort_snps.tranches \
    --recal-file output/cohort_snps.recal \
    -mode SNP

gatk --java-options "-Xmx7g" VariantFiltration \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/output.vqsr.vcf \
    -O output/output.vqsr.varfilter.vcf \
    --filter-name "Low_depth10" \
    --filter-expression "DP < 10"


 bcftools view -f 'PASS,.' -O vcf -o output/output.vqsr.varfilter.pass.vcf output/output.vqsr.varfilter.vcf

 bgzip -c output/output.vqsr.varfilter.pass.vcf > output/output.vqsr.varfilter.pass.vcf.gz
 tabix -p vcf output/output.vqsr.varfilter.pass.vcf.gz

 gatk VariantsToTable \
    -R reference/hg38/Homo_sapiens_assembly38.fasta \
    -V output/output.vqsr.varfilter.pass.vcf.gz \
    -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP \
    --show-filtered \
    -O output/output.vqsr.varfilter.pass.tsv
