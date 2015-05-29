#overall
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.var2vcf -r /home/lixia/ws/hg/hg19/human_g1k_v37.fasta -v example.var
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.var2vcf -r /home/lixia/ws/hg/hg19/human_g1k_v37.fasta -v example.var -f 1 -t bd
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.var2vcf -r /home/lixia/ws/hg/hg19/human_g1k_v37.fasta -v homo1.var -p homo1 -f 1 -t bd
