#overall
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -m test_fa.meta example.gap.bed simple.par example.fna
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -v test.lumpy example.gap.bed simple.par example.fna
PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -m test_fa.meta -v test.lumpy example.gap.bed simple.par example.fna
