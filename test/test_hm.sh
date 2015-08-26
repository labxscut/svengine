#in this example we generate one SV type of 1k size each in the region of 30M - 31M 
#how to run this with nohup?
nohup 2>&1 fasforge -o homo1.fas --outbam -n 1 -v homo1.var human.gap.bed simple.par human.fna >homo1.fas.out &
nohup 2>&1 mutforge -o homo1.mut --outbam -n 1 -v human.fna homo1.var human.gap.bed simple.par human.fna >homo1.mut.out &
