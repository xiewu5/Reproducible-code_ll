INDEX=/share/home/liuhy/01.genomics/35.Lates_calcarifer/01.onemore/02.Re-sequencing/03.ref/lc_genome.fa
cat config | while read id
do

arr=($id)
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
echo $sample $fq1 $fq2
bwa mem -t 100 -R "@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2 | samtools sort -@ 16 -o ../04.bam/$sample.bam -

done
