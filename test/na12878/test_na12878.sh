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

	### simulating using a reduced VAR file (.head.var) input one step into .bam file
	echo "### simulating NA12878 using a readuced size VAR file (.head.var)"
  PYTHONPATH=$SVENGINEPATH python -B -m mf.mutforge -x bam -n 32 -b 10000 -f 5000 -e 1000000 -v NA12878_1KG_SV.vcf.NA12878.head.var human_g1k_v37.fasta simple.par human_g1k_v37.fasta
  samtools index NA12878_1KG_SV.vcf.NA12878.head.lib1.bam

	### simulating using a full VAR file (.var) input fq files and map as a .bam file
  alias "swan_plot.R=/home/lixia/setup/swan/inst/swan_plot.R" #enable plotting
	echo "### simulating NA12878 using a VAR file (.var)"
  PYTHONPATH=$SVENGINEPATH python -B -m mf.mutforge -x fastq -n 32 -b 10000 -f 5000 -e 1000000 -v NA12878_1KG_SV.vcf.NA12878.var human_g1k_v37.fasta common.par human_g1k_v37.fasta
  bwa index human_g1k_v37.fasta
  bwa mem -t 32 human_g1k_v37.fasta NA12878_1KG_SV.vcf.NA12878.lib1.fq1 NA12878_1KG_SV.vcf.NA12878.lib1.fq2 | samtools view -@32 -O BAM -bh -o NA12878_1KG_SV.vcf.NA12878.lib1.bam -
  samtools sort -m 4G -@32 -o NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.bam
  samtools index NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam

  #simulate a null bam with 100x coverage; i.e. 1 billion 150*2 read pairs 
  wgsim -N 1000000000 -d 300 -s 5 -e 0.01 -r 0 -R 0 -X 0 -1 150 -2 150 human_g1k_v37.fasta NA12878_1KG_SV.null.NA12878.lib1.fq1 NA12878_1KG_SV.null.NA12878.lib1.fq2
  #nohup sh -c "wgsim -N 1000000000 -d 300 -s 5 -e 0.01 -r 0 -R 0 -X 0 -1 150 -2 150 human_g1k_v37.fasta NA12878_1KG_SV.null.NA12878.lib1.fq1 NA12878_1KG_SV.null.NA12878.lib1.fq2"
  bwa mem -t 32 human_g1k_v37.fasta NA12878_1KG_SV.null.NA12878.lib1.fq1 NA12878_1KG_SV.null.NA12878.lib1.fq2 | samtools view -@32 -O BAM -bh -o NA12878_1KG_SV.null.NA12878.lib1.bam -
  samtools sort -m 4G -@32 -o NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.bam
  samtools index NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam

	## output (vcf/null) files:
  # NA12878_1KG_SV.vcf.NA12878.fa  #fasta output, ...hap2.fa, ... if multiple haploids, contigs with spiked-in events of haploid 1
  # NA12878_1KG_SV.vcf.NA12878.fq1 #fastq output of read1, ...lib2.fq1, ... if have multiple libs sepcified in PAR (.par) file
  # NA12878_1KG_SV.vcf.NA12878.fq2 #fastq output of read2, ...lib2.fq2, ... if have multiple libs sepcified in PAR (.par) file
  # NA12878_1KG_SV.vcf.NA12878.lib1.bam #bam output, ...lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
  # NA12878_1KG_SV.vcf.NA12878.out.bed  #bed output for ground truth 
  # NA12878_1KG_SV.vcf.NA12878.out.var  #var output for ground truth 

  ## check events 
  head NA12878_1KG_SV.null.NA12878.out.bed >test_na12878.plot.bed
  swan_plot test_na12878.plot.bed NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam #looks OK
  swan_plot test_na12878.plot.bed NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam  #looks OK

  ## down-sampling to 90x, 75x, 50x, 25x, 10x
  samtools view -@32 -O BAM -bh -s 0.90 -o NA12878_1KG_SV.vcf.NA12878.lib1.90x.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.90x.sorted.bam 
  samtools view -@32 -O BAM -bh -s 0.75 -o NA12878_1KG_SV.vcf.NA12878.lib1.75x.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.75x.sorted.bam 
  samtools view -@32 -O BAM -bh -s 0.50 -o NA12878_1KG_SV.vcf.NA12878.lib1.50x.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.50x.sorted.bam
  samtools view -@32 -O BAM -bh -s 0.25 -o NA12878_1KG_SV.vcf.NA12878.lib1.25x.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.25x.sorted.bam
  samtools view -@32 -O BAM -bh -s 0.10 -o NA12878_1KG_SV.vcf.NA12878.lib1.10x.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.10x.sorted.bam
  
  ## down-sampling to 90x, 75x, 50x, 25x, 10x
  samtools view -@32 -O BAM -bh -s 0.90 -o NA12878_1KG_SV.null.NA12878.lib1.90x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.null.NA12878.lib1.90x.sorted.bam
  samtools view -@32 -O BAM -bh -s 0.75 -o NA12878_1KG_SV.null.NA12878.lib1.75x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.null.NA12878.lib1.75x.sorted.bam
  samtools view -@32 -O BAM -bh -s 0.50 -o NA12878_1KG_SV.null.NA12878.lib1.50x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.null.NA12878.lib1.50x.sorted.bam
  samtools view -@32 -O BAM -bh -s 0.25 -o NA12878_1KG_SV.null.NA12878.lib1.25x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.null.NA12878.lib1.25x.sorted.bam
  samtools view -@32 -O BAM -bh -s 0.10 -o NA12878_1KG_SV.null.NA12878.lib1.10x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && samtools index NA12878_1KG_SV.null.NA12878.lib1.10x.sorted.bam

  ## mix to 25%, 50% mixtures
  samtools view -H NA12878_1KG_SV.vcf.NA12878.lib1.100x.sorted.bam >inh.sam
  samtools merge -@32 -h inh.sam -O BAM NA12878_1KG_SV.vcf.NA12878.lib1.90%.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.90x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.10x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.90%.sorted.bam
  samtools merge -@32 -h inh.sam -O BAM NA12878_1KG_SV.vcf.NA12878.lib1.75%.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.75x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.25x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.75%.sorted.bam
  samtools merge -@32 -h inh.sam -O BAM NA12878_1KG_SV.vcf.NA12878.lib1.50%.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.50x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.50x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.50%.sorted.bam
  samtools merge -@32 -h inh.sam -O BAM NA12878_1KG_SV.vcf.NA12878.lib1.25%.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.25x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.75x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.25%.sorted.bam
  samtools merge -@32 -h inh.sam -O BAM NA12878_1KG_SV.vcf.NA12878.lib1.10%.sorted.bam NA12878_1KG_SV.vcf.NA12878.lib1.10x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.90x.sorted.bam && samtools index NA12878_1KG_SV.vcf.NA12878.lib1.10%.sorted.bam

  ## run lumpy, swan and theta

  ## compare sensitivity and specificity

  ## generate interval file for theta
  ## del
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.del.saf -o NA12878_1KG_SV.vcf.NA12878.homo.del.null.readcount.fct -T 32 NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && echo runtime=$SECONDS" >del.null.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.del.saf -o NA12878_1KG_SV.vcf.NA12878.homo.del.90%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.90%.sorted.bam && echo runtime=$SECONDS" >del.90%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.del.saf -o NA12878_1KG_SV.vcf.NA12878.homo.del.75%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.75%.sorted.bam && echo runtime=$SECONDS" >del.75%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.del.saf -o NA12878_1KG_SV.vcf.NA12878.homo.del.50%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.50%.sorted.bam && echo runtime=$SECONDS" >del.50%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.del.saf -o NA12878_1KG_SV.vcf.NA12878.homo.del.25%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.25%.sorted.bam && echo runtime=$SECONDS" >del.25%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.del.saf -o NA12878_1KG_SV.vcf.NA12878.homo.del.10%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.10%.sorted.bam && echo runtime=$SECONDS" >del.10%.count.log &
  ## dup
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.dup.saf -o NA12878_1KG_SV.vcf.NA12878.homo.dup.null.readcount.fct -T 32 NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && echo runtime=$SECONDS" >dup.null.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.dup.saf -o NA12878_1KG_SV.vcf.NA12878.homo.dup.90%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.90%.sorted.bam && echo runtime=$SECONDS" >dup.90%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.dup.saf -o NA12878_1KG_SV.vcf.NA12878.homo.dup.75%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.75%.sorted.bam && echo runtime=$SECONDS" >dup.75%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.dup.saf -o NA12878_1KG_SV.vcf.NA12878.homo.dup.50%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.50%.sorted.bam && echo runtime=$SECONDS" >dup.50%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.dup.saf -o NA12878_1KG_SV.vcf.NA12878.homo.dup.25%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.25%.sorted.bam && echo runtime=$SECONDS" >dup.25%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a NA12878_1KG_SV.vcf.NA12878.homo.dup.saf -o NA12878_1KG_SV.vcf.NA12878.homo.dup.10%.readcount.fct -T 32 NA12878_1KG_SV.vcf.NA12878.lib1.10%.sorted.bam && echo runtime=$SECONDS" >dup.10%.count.log &

  #separate ground truth events by four types: DEL, INS, DUP, INV
  grep -v "#" NA12878_1KG_SV.vcf.NA12878.groundtruth.bed | grep DEL >NA12878_1KG_SV.vcf.NA12878.groundtruth.del.bed #deletion
  grep -v "#" NA12878_1KG_SV.vcf.NA12878.groundtruth.bed | grep DUP >NA12878_1KG_SV.vcf.NA12878.groundtruth.dup.bed #duplication
  grep -v "#" NA12878_1KG_SV.vcf.NA12878.groundtruth.bed | grep INV >NA12878_1KG_SV.vcf.NA12878.groundtruth.inv.bed #inversion
  grep -v "#" NA12878_1KG_SV.vcf.NA12878.groundtruth.bed | grep "INS\|ALU\|SVA\|LINE" | grep -v "DEL" >NA12878_1KG_SV.vcf.NA12878.groundtruth.ins.bed #domestic and foreign insertion
  grep -v "#" NA12878_1KG_SV.vcf.NA12878.groundtruth.bed | grep "CNV" | grep human >>NA12878_1KG_SV.vcf.NA12878.groundtruth.dup.bed #CNV copy gain, go to DUP
  grep -v "#" NA12878_1KG_SV.vcf.NA12878.groundtruth.bed | grep "CNV" | grep -v human >>NA12878_1KG_SV.vcf.NA12878.groundtruth.del.bed #CNV copy loss, go to DEL
