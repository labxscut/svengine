#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -o homo1 -n 16 -v homo1.var human.fna human.gap.bed example.par human.fna
