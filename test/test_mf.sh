### A simple sanity checking script for mutforge ###
	# you can test even before installation
	# please run test_mf.sh /path/to/your/svengine/svengine
	if [ -z $1 ]; then 
		SVENGINEPATH=$HOME/setup/svengine/svengine/; 
	else
		SVENGINEPATH=$1; 
	fi
  
	## input files:
  # example.gap.bed #masked out regions in BED (.bed) format
  # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ... 
  # example.fna #nploid sequence fasta input, use comma to separate if multiple
  # example.fna #reference sequence fasta input
  # test_mf.meta #META file to specify event distribution
  # test_mf.var #VAR file to specify exact events

	## simulating using a META file (.meta) input
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 10000 -f 10000 -m test_mf.meta example.fna example.gap.bed simple.par example.fna
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge --outbam -n 16 -b 10000 -f 10000 -m test_mf.meta example.fna example.gap.bed simple.par example.fna

	## simulating using a VAR file (.var) input
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 10000 -f 10000 -v test_mf.var example.fna example.gap.bed simple.par example.fna
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge --outbam -n 16 -b 10000 -f 10000 -v test_mf.var example.fna example.gap.bed simple.par example.fna

	## simulating using both a VAR file (.var) and a META file (.meta) input (experimental)
  #PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 100 -f 100 -m test_mf.meta -v test_mf.var example.fna example.gap.bed simple.par example.fna

	## output files:
  # test_mf.lib1.bam #bam output, test_mf.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.out.bed #bed output for ground truth 
  # test_mf.out.var #var output for ground truth 
