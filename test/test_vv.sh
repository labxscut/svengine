#testing components
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -m DEL.meta example.gap.bed simple.par example.fna
