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
  # test_mf.meta #META file to specify event distribution
  # test_mf.var #VAR file to specify exact events

  # -f 10000: plansize 10000
  # NOTE: choose your plansize wisely to provide best randomization, plansize * number_of_variants * [10 ~ 100] = genomesize 
  # -b 100000: burnin 100000
  # NOTE: choose your burnin wisely to avoid telomere N regions, N ~= 1000000 for human genome
  # -e 100000: trunksize 100000
  # NOTE: choose your trunksize wisely to save running time such that trunksize*nprocs~=genomesize

	### testing using a META file (.meta) input or a VAR file (.var) input
	echo "### testing using a META file (.meta) input"
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
	echo "### testing using a VAR file (.var) input"
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 1 -b 100000 -f 10000 -e 100000 -v test_mf.var example.fna simple.par example.fna

  ### testing with multiple (16) process
  echo "### testing with multiple (16) process"
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna

  ### testing output fasta, fastq or bam
  echo "### testing output fasta"
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -x fasta -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
  echo "### testing output fastq"
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -x fastq -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
  echo "### testing output bam"
  PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -x bam -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
  samtools index test_mf.lib1.bam

	## output files:
  # test_mf.hap1.fa  #fasta output, ...hap2.fa, ... if multiple haploids, contigs with spiked-in events of haploid 1
  # test_mf.lib1.fq1 #fastq output of read1, ...lib2.fq1, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.lib1.fq2 #fastq output of read2, ...lib2.fq2, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.lib1.bam #bam output, ...lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.out.bed  #bed output for ground truth 
  # test_mf.out.var  #var output for ground truth 
