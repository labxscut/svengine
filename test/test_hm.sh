#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
#how to run this with nohup?
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -o homo1.fas -n 1 -v homo1.var human.gap.bed example.par human.fna >homo1.fas.out 2>&1 &
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge -o homo1.mut -n 1 -v homo1.var human.fna human.gap.bed example.par human.fna >homo1.mut.out 2>&1 &
