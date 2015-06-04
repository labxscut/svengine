### A simple sanity checking script for fasforge ###
  
	## input files:
  # example.gap.bed #masked out regions in BED (.bed) format
  # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ... 
  # example.fna #fasta input 
  # test_fa.meta #META file to specify event distribution
  # test_fa.var #VAR file to specify exact events

	## simulating using a META file (.meta) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -m test_fa.meta example.gap.bed simple.par example.fna

	## simulating using a VAR file (.var) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -v test_fa.var example.gap.bed simple.par example.fna

	## simulating using both a VAR file (.var) and a META file (.meta) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.fasforge -m test_fa.meta -v test_fa.var example.gap.bed simple.par example.fna

	## output files:
  # test_fa.lib1.bam #bam output, test_fa.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # test_fa.out.bed #bed output for ground truth 
  # test_fa.out.var #var output for ground truth 
  # test_fa.fnv #fasta output after ground truth spiked-in 
