# sort the bam by read name
singularity exec ~/biosoft/Assembly202306.sif samtools sort -n -@ 12 -o aligned.sort_name.bam  sample.bwa_mem.bam

# matlock bam2juicer

singularity exec ~/biosoft/Assembly202306.sif  matlock bam2 juicer aligned.sort_name.bam   merged_nodups.txt

# agp to assembly
singularity exec ~/biosoft/Assembly202306.sif  agp2assembly.py  groups.agp    groups.assembly

# juicetools
singularity exec ~/biosoft/Assembly202306.sif \
  bash /opt/3d-dna/visualize/run-assembly-visualizer.sh  -q 1 -p true groups.assembly  merged_nodups.txt

