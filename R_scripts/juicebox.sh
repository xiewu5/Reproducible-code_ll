#!/bin/bash

ln -s /share/home/liuhy/01.genomics/01.Channa_argus/71.xy/my_T2T_project/03_hic_scaffolding/ca_genomics.fa .
samtools faidx ca_genomics.fa
/share/home/liuhy/biosoft/HapHiC/scripts/../utils/juicer pre -a -q 1 -o out_JBAT /share/home/liuhy/01.genomics/01.Channa_argus/71.xy/my_T2T_project/03_hic_scaffolding/HiC.filtered.bam scaffolds.raw.agp ca_genomics.fa.fai >out_JBAT.log 2>&1
(java -Djava.awt.headless=true -jar -Xmx320G /share/home/liuhy/biosoft/HapHiC/scripts/../utils/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
