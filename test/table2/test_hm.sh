### A simple sanity checking script for mutforge with human reference hg19 ###
	# you can test even before installation
	# please run test_hm.sh /path/to/your/svengine/svengine
	if [ -z $1 ]; then 
		SVENGINEPATH=$HOME/setup/svengine/svengine/; 
	else
		SVENGINEPATH=$1; 
	fi
  
	## input files:
  # humanbed #masked out regions in BED (.bed) format
  # test_hm.par #PAR file to specify library type: insert distribution, coverage, read length, ... 
  # example.fna #nploid sequence fasta input, use comma to separate if multiple
  # test_hm.meta #META file to specify event distribution
  # test_hm.var #VAR file to specify exact events

	## simulating upto FASTA output
  nohup sh -c 'export PYTHONPATH=$SVENGINEPATH; python -m mf.mutforge -b 0 -e 1000000 -x fasta -o "test_hm.nproc=1.output=fasta" -n 1 -v test_hm.var human_g1k_v37.fasta test_hm.par human_g1k_v37.fasta >"test_hm.nproc=1.output=fasta.log" 2>&1' &
  nohup sh -c 'export PYTHONPATH=$SVENGINEPATH; python -m mf.mutforge -b 0 -e 1000000 -x fasta -o "test_hm.nproc=64.output=fasta" -n 64 -v test_hm.var human_g1k_v37.fasta test_hm.par human_g1k_v37.fasta >"test_hm.nproc=64.output=fasta.log" 2>&1' &

	## simulating upto FASTQ output
  nohup sh -c 'export PYTHONPATH=$SVENGINEPATH; python -m mf.mutforge -b 0 -e 1000000 -x fastq -o "test_hm.nproc=1.output=fastq" -n 1 -v test_hm.var human_g1k_v37.fasta test_hm.par human_g1k_v37.fasta >"test_hm.nproc=1.output=fastq.log" 2>&1' &
  nohup sh -c 'export PYTHONPATH=$SVENGINEPATH; python -m mf.mutforge -b 0 -e 1000000 -x fastq -o "test_hm.nproc=64.output=fastq" -n 64 -v test_hm.var human_g1k_v37.fasta test_hm.par human_g1k_v37.fasta >"test_hm.nproc=64.output=fastq.log" 2>&1' &

  ## simultating upto BAM output
  nohup sh -c 'export PYTHONPATH=$SVENGINEPATH; python -m mf.mutforge -b 0 -e 1000000 -x bam -o "test_hm.nproc=1.output=bam" -n 1 -v test_hm.var human_g1k_v37.fasta test_hm.par human_g1k_v37.fasta >"test_hm.nproc=1.output=bam.log" 2>&1' &
  nohup sh -c 'export PYTHONPATH=$SVENGINEPATH; python -m mf.mutforge -b 0 -e 1000000 -x bam -o "test_hm.nproc=64.output=bam" -n 64 -v test_hm.var human_g1k_v37.fasta test_hm.par human_g1k_v37.fasta >"test_hm.nproc=64.output=bam.log" 2>&1' &
