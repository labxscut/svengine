#this is single threaded tests on human_g1k_v37.fasta
nohup 2>&1 fasforge -o homo1.fas --outbam -n 1 -v homo1.var human.gap.bed simple.par human_g1k_v37.fasta >homo1.fas.out &
nohup 2>&1 mutforge -o homo1.mut --outbam -n 1 -v homo1.var human_g1k_v37.fasta human.gap.bed simple.par human_g1k_v37.fasta >homo1.mut.out &
