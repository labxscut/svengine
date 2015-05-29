#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -o control -n 16 -m control.meta example.fna example.gap.bed example.par example.fna
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -o tumor -n 16 -m tumor.meta example.gap.bed example.par example.fna
