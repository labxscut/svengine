#overall
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -m example.meta example.gap.bed simple.par example.fna
PYTHONPATH=$HOME/setup/svengine/mutforge/ python -m mf.fasforge -v test.lumpy example.gap.bed simple.par example.fna
