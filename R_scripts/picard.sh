ls *bam | cut -d"." -f 1 | sort -u | while read id; do

picard -Xmx320g MarkDuplicates I=${id}.bam  O=${id}.mkdup.bam  CREATE_INDEX=true M=${id}.marked_dup_metrics.txt

done
