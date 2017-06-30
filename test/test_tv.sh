### A simple sanity checking script for tree2var ###
	# you can test even before installation
	# please run test_tv.sh /path/to/your/svengine/svengine
	if [ -z $1 ]; then 
		SVENGINEPATH=$HOME/setup/svengine/svengine/; 
	else
		SVENGINEPATH=$1; 
	fi
  
  ## Step 1 - prepare tree-based VAR file

	## input files:
  # egTree.newick #Newick formatted phylogenetic tree input
  # egTree.var    #VAR file to specify exact events, whose frequency will be changed based on tree

	## preparing a new VAR file based on egTree frequencies input
  PYTHONPATH=$SVENGINEPATH python -m mf.tree2var -p egTree.var egTree.newick

	## output files:
  # egTree.tree.var #VAR output to be used as input for simultation by mutforge
  # egTree.tree.plot #a text plot of the assigned phylogenetic tree

  ## Step 2 - simulate tree-based VAR file

	## input files:
  # egTree.tree.var #VAR output to be used as input for simultation by mutforge
  # example.gap.bed #masked out regions in BED (.bed) format
  # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ...
  # example.fna #nploid sequence fasta input, use comma to separate if multiple

	## simulate clonal evolution using tree-based VAR file by mutforge
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge --outbam -n 16 -b 10000 -f 10000 -e 10000 -v egTree.tree.var -s example.gap.bed example.fna simple.par example.fna

  ## output files:
  # egTree.tree.lib1.bam #bam output, test_mf.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # egTree.tree.lib1.fq1/2 #fastq output, test_mf.lib2.fq1/2, ... if have multiple libs sepcified in PAR (.par) file
  # egTree.tree.out.bed #bed output for ground truth
