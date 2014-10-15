#overall
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -m test_homo.meta -v test_homo.lumpy example.bed simple.par example.fna
#test deletion
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -v fadel.lumpy example.bed simple.par example.fna
#test duplication
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -v fadup.lumpy example.bed simple.par example.fna
#test translocation
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -v fatra.lumpy example.bed simple.par example.fna
#test inversion
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -v fainv.lumpy example.bed simple.par example.fna
#test finsertion
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -m fafins.meta example.bed simple.par example.fna
#test dinsertion
PYTHONPATH=$HOME/setup/mutforge/mutforge/ python -m mf.fasforge -m fadins.meta example.bed simple.par example.fna
