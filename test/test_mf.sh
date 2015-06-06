### A simple sanity checking script for mutforge ###
  
	## input files:
  # example.gap.bed #masked out regions in BED (.bed) format
  # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ... 
  # example.fna #nploid sequence fasta input, use comma to separate if multiple
  # example.fna #reference sequence fasta input
  # test_mf.meta #META file to specify event distribution
  # test_mf.var #VAR file to specify exact events

	## simulating using a META file (.meta) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge -n 16 -f 100 -m test_mf.meta example.fna example.gap.bed simple.par example.fna

	## simulating using a VAR file (.var) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge -n 16 -f 100 -v test_mf.var example.fna example.gap.bed simple.par example.fna

	## simulating using both a VAR file (.var) and a META file (.meta) input
  PYTHONPATH=$HOME/setup/svengine/svengine/ python -m mf.mutforge -n 16 -f 100 -m test_mf.meta -v test_mf.var example.fna example.gap.bed simple.par example.fna

	## output files:
  # test_mf.lib1.bam #bam output, test_mf.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.out.bed #bed output for ground truth 
  # test_mf.out.var #var output for ground truth 
