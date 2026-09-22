#!/bin/bash
# run_pipeline.sh
# Full command sequence for the NA12878 chromosome 20 variant-calling pipeline.
# Run stages individually and verify output at each step — see README.md for
# the reasoning behind each stage. Byte offsets in Stage 1 are specific to
# this run (calculated from the reference's .fai index for chr20).
set -e  # stop immediately if any command fails.

# ---------------------------------------------------------------------------
# Stage 0: Environments
# ---------------------------------------------------------------------------
conda create -n chr20_variant_calling python=3.10 sra-tools samtools bwa fastqc fastp gatk4 -c bioconda -c conda-forge -y
conda create -n happy_env -c bioconda -c conda-forge hap.py rtg-tools -y
conda create -n snpeff_env -c bioconda -c conda-forge snpeff openjdk=21 -y

conda activate chr20_variant_calling

# ---------------------------------------------------------------------------
# Stage 1: Extract chr20 reference via byte-range request
# ---------------------------------------------------------------------------
curl -O "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
grep "^chr20" GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
# chr20  64444167  2751788762  70  71  (length, start_offset, bases/line, bytes/line)

curl -r 2751788762-2817153559 \
  "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
  -o chr20_body.fa

echo ">chr20" > chr20_reference.fa
cat chr20_body.fa >> chr20_reference.fa
samtools faidx chr20_reference.fa

# ---------------------------------------------------------------------------
# Stage 2: Extract chr20 reads from the remote CRAM
# ---------------------------------------------------------------------------
samtools view -b -T chr20_reference.fa \
  -o na12878_chr20.bam \
  "ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram" chr20

samtools flagstat na12878_chr20.bam

# ---------------------------------------------------------------------------
# Stage 3: Convert back to raw paired FASTQ
# ---------------------------------------------------------------------------
samtools sort -n na12878_chr20.bam -o na12878_chr20.namesorted.bam
samtools fastq \
  -1 na12878_chr20_R1.fastq -2 na12878_chr20_R2.fastq \
  -0 /dev/null -s /dev/null \
  na12878_chr20.namesorted.bam

# ---------------------------------------------------------------------------
# Stage 4: QC and trimming
# ---------------------------------------------------------------------------
fastqc na12878_chr20_R1.fastq na12878_chr20_R2.fastq

fastp \
  -i na12878_chr20_R1.fastq -I na12878_chr20_R2.fastq \
  -o na12878_chr20_R1.trimmed.fastq -O na12878_chr20_R2.trimmed.fastq \
  -h fastp_report.html -j fastp_report.json

fastqc na12878_chr20_R1.trimmed.fastq na12878_chr20_R2.trimmed.fastq

# ---------------------------------------------------------------------------
# Stage 5: Alignment
# ---------------------------------------------------------------------------
bwa index chr20_reference.fa
bwa mem -t 4 chr20_reference.fa \
  na12878_chr20_R1.trimmed.fastq na12878_chr20_R2.trimmed.fastq \
  > na12878_chr20.aligned.sam

samtools view -b na12878_chr20.aligned.sam > na12878_chr20.aligned.bam

# ---------------------------------------------------------------------------
# Stage 6: Duplicate marking (requires name-sort + fixmate first)
# ---------------------------------------------------------------------------
samtools sort -n na12878_chr20.aligned.bam -o na12878_chr20.namesorted_for_fixmate.bam
samtools fixmate -m na12878_chr20.namesorted_for_fixmate.bam na12878_chr20.fixmate.bam
samtools sort na12878_chr20.fixmate.bam -o na12878_chr20.sorted.bam
samtools index na12878_chr20.sorted.bam
samtools markdup na12878_chr20.sorted.bam na12878_chr20.dedup.bam

samtools flagstat na12878_chr20.dedup.bam

# ---------------------------------------------------------------------------
# Stage 7: Read groups + sequence dictionary
# ---------------------------------------------------------------------------
samtools addreplacerg \
  -r "@RG\tID:na12878_chr20\tSM:NA12878\tLB:lib1\tPL:ILLUMINA" \
  -o na12878_chr20.rg.bam na12878_chr20.dedup.bam
samtools index na12878_chr20.rg.bam
samtools dict chr20_reference.fa > chr20_reference.dict

# ---------------------------------------------------------------------------
# Stage 8: Variant calling (GATK HaplotypeCaller)
# ---------------------------------------------------------------------------
gatk --java-options "-Xmx2g" HaplotypeCaller \
  -R chr20_reference.fa -I na12878_chr20.rg.bam \
  -O na12878_chr20.raw.vcf.gz

# ---------------------------------------------------------------------------
# Stage 9: Split by type and apply hard filters
# ---------------------------------------------------------------------------
gatk SelectVariants -R chr20_reference.fa -V na12878_chr20.raw.vcf.gz \
  --select-type-to-include SNP -O na12878_chr20.snps.vcf.gz
gatk SelectVariants -R chr20_reference.fa -V na12878_chr20.raw.vcf.gz \
  --select-type-to-include INDEL -O na12878_chr20.indels.vcf.gz

gatk VariantFiltration -R chr20_reference.fa -V na12878_chr20.snps.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O na12878_chr20.snps.filtered.vcf.gz

gatk VariantFiltration -R chr20_reference.fa -V na12878_chr20.indels.vcf.gz \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O na12878_chr20.indels.filtered.vcf.gz

gatk MergeVcfs \
  -I na12878_chr20.snps.filtered.vcf.gz -I na12878_chr20.indels.filtered.vcf.gz \
  -O na12878_chr20.filtered.vcf.gz

# ---------------------------------------------------------------------------
# Stage 10: Truth set preparation
# ---------------------------------------------------------------------------
curl -O "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
curl -O "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
curl -O "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"

bcftools view -r chr20 HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O z -o truth_chr20.vcf.gz
tabix -p vcf truth_chr20.vcf.gz
awk '$1 == "chr20"' HG001_GRCh38_1_22_v4.2.1_benchmark.bed > confident_chr20.bed

# ---------------------------------------------------------------------------
# Stage 11: Benchmarking (separate environment — hap.py needs Python 2)
# ---------------------------------------------------------------------------
# conda activate happy_env
export HGREF=chr20_reference.fa
hap.py truth_chr20.vcf.gz na12878_chr20.filtered.vcf.gz \
  -f confident_chr20.bed -r chr20_reference.fa \
  --engine=vcfeval -o benchmark_results

# ---------------------------------------------------------------------------
# Stage 12: Annotation (separate environment — SnpEff needs Java 21)
# ---------------------------------------------------------------------------
# conda activate snpeff_env
snpEff download GRCh38.p14
export _JAVA_OPTIONS="-Xmx3g"
snpEff GRCh38.p14 na12878_chr20.filtered.vcf.gz > na12878_chr20.annotated.vcf

grep -v "^#" na12878_chr20.annotated.vcf | grep "HIGH" > high_impact_variants.txt

# ---------------------------------------------------------------------------
# Stage 13: Result figures
# ---------------------------------------------------------------------------
# conda activate chr20_variant_calling
python3 scripts/plot_results.py