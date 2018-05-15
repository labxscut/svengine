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
  echo "### simulating theta example using Exmaple.intervals derived VAR file"
  nohup sh -c "PYTHONPATH=$SVENGINEPATH python -B -m mf.mutforge -x bam -n 64 -b 10000 -f 5000 -e 10000 -v Example.intervals.var human_g1k_v37.fasta common.par human_g1k_v37.fasta && samtools index Example.intervals.lib1.bam && echo $SECONDS" >Example.intervals.log

  ## down-sampling to 90x, 75x, 50x, 25x, 10x and generating mixtures
  samtools view -H Example.intervals.lib1.100x.sorted.bam >inh.sam
  nohup sh -c "samtools view -@64 -O BAM -bh -s 0.90 -o Example.intervals.lib1.090x.sorted.bam Example.intervals.lib1.100x.sorted.bam && samtools index Example.intervals.lib1.090x.sorted.bam && samtools merge -@64 -h inh.sam -O BAM Example.intervals.lib1.90%.sorted.bam Example.intervals.lib1.090x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.10x.sorted.bam && samtools index Example.intervals.lib1.90%.sorted.bam" >theta.90%.out
  nohup sh -c "samtools view -@64 -O BAM -bh -s 0.75 -o Example.intervals.lib1.075x.sorted.bam Example.intervals.lib1.100x.sorted.bam && samtools index Example.intervals.lib1.075x.sorted.bam && samtools merge -@64 -h inh.sam -O BAM Example.intervals.lib1.75%.sorted.bam Example.intervals.lib1.075x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.25x.sorted.bam && samtools index Example.intervals.lib1.75%.sorted.bam" >theta.75%.out
  nohup sh -c "samtools view -@64 -O BAM -bh -s 0.50 -o Example.intervals.lib1.050x.sorted.bam Example.intervals.lib1.100x.sorted.bam && samtools index Example.intervals.lib1.050x.sorted.bam && samtools merge -@64 -h inh.sam -O BAM Example.intervals.lib1.50%.sorted.bam Example.intervals.lib1.050x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.50x.sorted.bam && samtools index Example.intervals.lib1.50%.sorted.bam" >theta.50%.out
  nohup sh -c "samtools view -@64 -O BAM -bh -s 0.25 -o Example.intervals.lib1.025x.sorted.bam Example.intervals.lib1.100x.sorted.bam && samtools index Example.intervals.lib1.025x.sorted.bam && samtools merge -@64 -h inh.sam -O BAM Example.intervals.lib1.25%.sorted.bam Example.intervals.lib1.025x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.75x.sorted.bam && samtools index Example.intervals.lib1.25%.sorted.bam" >theta.25%.out
  nohup sh -c "samtools view -@64 -O BAM -bh -s 0.10 -o Example.intervals.lib1.010x.sorted.bam Example.intervals.lib1.100x.sorted.bam && samtools index Example.intervals.lib1.010x.sorted.bam && samtools merge -@64 -h inh.sam -O BAM Example.intervals.lib1.10%.sorted.bam Example.intervals.lib1.010x.sorted.bam NA12878_1KG_SV.null.NA12878.lib1.90x.sorted.bam && samtools index Example.intervals.lib1.10%.sorted.bam" >theta.10%.out

  ## counting features

  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a Example.intervals.saf -o Example.intervals.lib1.10%.sorted.bam.readcount.fct -T 32 Example.intervals.lib1.10%.sorted.bam && echo runtime=$SECONDS" >theta.10%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a Example.intervals.saf -o Example.intervals.lib1.25%.sorted.bam.readcount.fct -T 32 Example.intervals.lib1.25%.sorted.bam && echo runtime=$SECONDS" >theta.25%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a Example.intervals.saf -o Example.intervals.lib1.50%.sorted.bam.readcount.fct -T 32 Example.intervals.lib1.50%.sorted.bam && echo runtime=$SECONDS" >theta.50%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a Example.intervals.saf -o Example.intervals.lib1.75%.sorted.bam.readcount.fct -T 32 Example.intervals.lib1.75%.sorted.bam && echo runtime=$SECONDS" >theta.75%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a Example.intervals.saf -o Example.intervals.lib1.90%.sorted.bam.readcount.fct -T 32 Example.intervals.lib1.90%.sorted.bam && echo runtime=$SECONDS" >theta.90%.count.log &
  nohup sh -c "$HOME/setup/subread-1.6.1-source/bin/featureCounts -F SAF -a Example.intervals.saf -o Example.intervals.lib1.00%.sorted.bam.readcount.fct -T 32 NA12878_1KG_SV.null.NA12878.lib1.100x.sorted.bam && echo runtime=$SECONDS" >theta.00%.count.log &
