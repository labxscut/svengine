### A script for using mutforge to embed 20863 NA12878 events ###
	if [ -z $1 ]; then 
		SVENGINEPATH=$HOME/setup/svengine/svengine/; 
	else
		SVENGINEPATH=$1; 
	fi
  
	## input files:
  # example.gap.bed #masked out regions in BED (.bed) format
  # NA12878_1KG_SV.vcf.NA12878.var #VAR file to specify exact events
  # human_g1k_v37.fasta #template hg37 sequence fasta input
  # common.par #PAR file to specify common single lib design: insert distribution, coverage, read length, ... 

  # -f 5000: plansize 5000
  # NOTE: size of basic parallelized block
  # -b 100000: burnin region without simulation
  # NOTE: choose your burnin wisely to avoid telomere N regions, N ~= 1000000 for human genome
  # -e 1000000: trunksize 1000000
  # NOTE: choose your trunksize wisely to save running time such that trunksize*nprocs~=genomesize

	### simulating using a VAR file (.var) input
	echo "### simulating NA12879 using a VAR file (.var)"
  PYTHONPATH=$SVENGINEPATH python -B -m mf.mutforge -n 16 -b 10000 -f 5000 -e 1000000 -v NA12878_1KG_SV.vcf.NA12878.var human_g1k_v37.fasta common.par human_g1k_v37.fasta
  PYTHONPATH=$SVENGINEPATH python -B -m mf.mutforge -x bam -n 16 -b 10000 -f 5000 -e 1000000 -v NA12878_1KG_SV.vcf.NA12878.head.var human_g1k_v37.fasta simple.par human_g1k_v37.fasta
  samtools index test_mf.lib1.bam

	## output files:
  # test_mf.hap1.fa  #fasta output, ...hap2.fa, ... if multiple haploids, contigs with spiked-in events of haploid 1
  # test_mf.lib1.fq1 #fastq output of read1, ...lib2.fq1, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.lib1.fq2 #fastq output of read2, ...lib2.fq2, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.lib1.bam #bam output, ...lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # test_mf.out.bed  #bed output for ground truth 
  # test_mf.out.var  #var output for ground truth 
