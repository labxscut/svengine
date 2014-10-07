#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
#PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge -m example.meta example.fasta example.bed example.par example.fasta
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge --nprocs 1  -m deletion.meta example.fna example.bed deletion.par example.fna
