#!/bin/bash

ref=../03.move_redundancy/genome.final.fasta
R1=../00.data/hic_R1.clean.fq.gz
R2=../00.data/hic_R2.clean.fq.gz
threads=100

group_count=24
enzyme=GATC  # HindIII: AAGCTT; MboI: GATC

### index reference genome
ln -s ${ref}  ./seq.HiCcorrected.fasta 
bwa index seq.HiCcorrected.fasta
samtools faidx seq.HiCcorrected.fasta

bwa mem -SP5M -t $threads seq.HiCcorrected.fasta $R1 $R2 \
  | samtools  view -hF 256 - \
  | samtools  sort -@ $threads -o sample.bwa_mem.bam -T tmp.ali


### filter bam
samtools view -bq 40 sample.bwa_mem.bam | \
  samtools view  -bt seq.HiCcorrected.fasta.fai > sample.unique.bam
PreprocessSAMs.pl sample.unique.bam seq.HiCcorrected.fasta $enzyme

### partition
ALLHiC_partition \
  -r seq.HiCcorrected.fasta \
  -e $enzyme \
  -k $group_count \
  -b sample.unique.REduced.paired_only.bam

### optimize
rm cmd.list
for((K=1;K<=$group_count;K++));do echo "allhic optimize \
  sample.unique.REduced.paired_only.counts_${enzyme}.${group_count}g${K}.txt \
  sample.unique.REduced.paired_only.clm" >> cmd.list;done

ParaFly -c cmd.list -CPU $threads

### build
ALLHiC_build  seq.HiCcorrected.fasta

### plot
samtools faidx groups.asm.fasta
cut -f1,2 groups.asm.fasta.fai|grep sample > chrn.list
ALLHiC_plot sample.bwa_mem.bam groups.agp chrn.list 500k pdf


