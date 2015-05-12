#overall
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -m test.meta -v test.lumpy example.gap.bed simple.par example.fna
#test deletion
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -v fadel.lumpy example.gap.bed simple.par example.fna
#test duplication
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -v fadup.lumpy example.gap.bed simple.par example.fna
#test translocation
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -v fatra.lumpy example.gap.bed simple.par example.fna
#test inversion
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -v fainv.lumpy example.gap.bed simple.par example.fna
#test finsertion
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -m fafins.meta example.gap.bed simple.par example.fna
#test dinsertion
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -m fadins.meta example.gap.bed simple.par example.fna
#test nprocs
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -n 2 -q 1,0.5 -m fadins.meta example.gap.bed simple.par example.fna >simple.log
#test scaling nprocs multilib multifreq
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -n 10 -q 1,0.5 -m fadins.meta example.gap.bed first.par example.fna >scale.log
