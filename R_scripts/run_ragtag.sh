ln -s  ../02.RMcontaminate/genome.final.fasta 
singularity exec ~/biosoft/Assembly202306.sif  ragtag.py  scaffold \
  -u -t 100 ../data/ref.fasta  ./genome.final.fasta
