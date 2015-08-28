#!/bin/bash

#A simple script to test if embedding SV regions are good
#./test_fa.sh #generating bamfiles
#./test_mf.sh #generating bamfiles
samtools sort test_fa.lib1.bam test_fa.st
samtools sort test_mf.lib1.bam test_mf.st
samtools index test_fa.st.bam
samtools index test_mf.st.bam
#plotting out embedded regions, require installation of SWAN@http://bitbucket.org/charade/swan
swan_plot.R -f bed -o test_fa.out test_fa.out.bed test_fa.st.bam
swan_plot.R -f bed -o test_mf.out test_mf.out.bed test_mf.st.bam

#### A simple sanity checking script for fasforge ### 
#  ## input files:
#  # example.gap.bed #masked out regions in BED (.bed) format
#  # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ...
#  # example.fna #fasta input
#  # test_fa.meta #META file to specify event distribution
#  # test_fa.var #VAR file to specify exact events
#
#  ## simulating using a META file (.meta) input
#  fasforge -n 16 -b 100 -f 100 -m test_fa.meta example.gap.bed simple.par example.fna
#
#  ## simulating using a VAR file (.var) input
#  fasforge -n 16 -b 100 -f 100 -v test_fa.var example.gap.bed simple.par example.fna
#
#  ## simulating using both a VAR file (.var) and a META file (.meta) input (experimental)
#  #PYTHONPATH=$SVENGINEPATH python -m mf.fasforge -n 16 -b 100 -f 100 -m test_fa.meta -v test_fa.var example.gap.bed simple.par example.fna
#
#  ## output files:
#  # test_fa.lib1.bam #bam output, test_fa.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
#  # test_fa.out.bed #bed output for ground truth
#  # test_fa.out.var #var output for ground truth
#  # test_fa.fnv #fasta output after ground truth spiked-in
#
#### A simple sanity checking script for mutforge ###
#
#  ## input files:
#  # example.gap.bed #masked out regions in BED (.bed) format
#  # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ...
#  # example.fna #nploid sequence fasta input, use comma to separate if multiple
#  # example.fna #reference sequence fasta input
#  # test_mf.meta #META file to specify event distribution
#  # test_mf.var #VAR file to specify exact events
#
#  ## simulating using a META file (.meta) input
#  mutforge -n 16 -b 100 -f 100 -m test_mf.meta example.fna example.gap.bed simple.par example.fna
#
#  ## simulating using a VAR file (.var) input
#  mutforge -n 16 -b 100 -f 100 -v test_mf.var example.fna example.gap.bed simple.par example.fna
#
#  ## simulating using both a VAR file (.var) and a META file (.meta) input (experimental)
#  #PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 100 -f 100 -m test_mf.meta -v test_mf.var example.fna example.gap.bed simple.par example.fna
#
#  ## output files:
#  # test_mf.lib1.bam #bam output, test_mf.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
#  # test_mf.out.bed #bed output for ground truth
#  # test_mf.out.var #var output for ground truth
