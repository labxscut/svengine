### A simple sanity checking script for fasforge ###
  
	## simulating using a META file (.meta) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -m test_fa.meta example.gap.bed simple.par example.fna

	## simulating using a VAR file (.var) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -v test_fa.var example.gap.bed simple.par example.fna

	## simulating using both a VAR file (.var) and a META file (.meta) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -m test_fa.meta -v test_fa.var example.gap.bed simple.par example.fna
