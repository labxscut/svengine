#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
#python -m ../mutforge.mf.mutforge.py -m example.meta example.fasta example.bed example.par
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.mutforge -m example.meta example.fasta example.bed example.par example.fasta
