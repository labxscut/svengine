.. |TREE| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/TREE.png
   :alt: TREE.png
   :width: 800px

.. contents:: Navigator

Examples
========

To easy the initial learning, we include multiple BASH scripts of SVEngine .
The scripts represent common analysis tasks with prepared input files. 

To run following examples:

  - Download and unzip SVEngine package or clone the source repository to $SVENGINEPATH.
  - Export the $SVENGINEPATH to your SHELL: 'export SVENGINEPATH=$SVENGINEPATH'
  - Change directory into $SVENGINEPATH/test and run the test_dep.sh to make sure dependencies OK
  - run test_mf.sh and test_tv.sh scripts for example tasks
  - required input and output files were explained inline withe scripts and in this wiki
  - run clean.sh to clear outputs from previous runs

An example SVEngine mutforge run 
--------------------------------------------------

To run this example, run the 'test_mf.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir: 

$SVENGINEPATH/test/test_mf.sh

  .. code:: bash

    #!/bin/bash
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
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
    	echo "### testing using a VAR file (.var) input"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -n 1 -b 100000 -f 10000 -e 100000 -v test_mf.var example.fna simple.par example.fna
    
      ### testing with multiple (16) process
      echo "### testing with multiple (16) process"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -n 16 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
    
      ### testing output fasta, fastq or bam
      echo "### testing output fasta"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -x fasta -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
      echo "### testing output fastq"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -x fastq -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
      echo "### testing output bam"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -x bam -n 1 -b 100000 -f 10000 -e 100000 -m test_mf.meta example.fna simple.par example.fna
      samtools index test_mf.lib1.bam
    
    	## output files:
      # test_mf.hap1.fa  #fasta output, ...hap2.fa, ... if multiple haploids, contigs with spiked-in events of haploid 1
      # test_mf.lib1.fq1 #fastq output of read1, ...lib2.fq1, ... if have multiple libs sepcified in PAR (.par) file
      # test_mf.lib1.fq2 #fastq output of read2, ...lib2.fq2, ... if have multiple libs sepcified in PAR (.par) file
      # test_mf.lib1.bam #bam output, ...lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
      # test_mf.out.bed  #bed output for ground truth 
      # test_mf.out.var  #var output for ground truth 


An example SVEngine tree2var run 
--------------------------------------------------

To run this example, run the 'test_tv.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir: 

$SVENGINEPATH/test/test_tv.sh

  .. code:: bash

    #!/bin/bash
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
      # egTree.newick #Newick formatted phylogenetic tree input with special restrictions below:
                      #The input must represent a binary tree with each internal node named by 
                      #one variant's meta id (MID), 
                      #and each leaf node named by one of cell subpopulations.
                      #the first of the only two children of an internal node always represent
                      #the cell populaiton HAS NOT gained the mutation denoted by the internal node.  
                      #the second of the only two children of an internal node always represent
                      #the cell populaiton HAS gained the mutation denoted by the internal node.  
      # egTree.var    #VAR file to specify exact mutations, whose frequency will be overwritten 
                      #based on the phylogeny information provided in egTree.newick
                      #The variants meta id (MID) in this list must match to those
                      #in egTree.newick.
                      #Optionally, this VAR file can be generated by fasforge or mutforge with --layout
    
    	### preparing a new VAR file based on egTree frequencies input
      echo "### preparing a new VAR file based on egTree frequencies input"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.tree2var -p egTree.var egTree.newick
    
    	## output files:
      # egTree.tree.var #VAR output to be used as input for simultation by mutforge
      # egTree.tree.plot #a text plot of the assigned phylogenetic tree
      # egTree.tree.json #a json output of variant frquency calculation of the assigned phylogenetic tree
    
      ## Step 2 - simulate tree-based VAR file
    
    	## input files:
      # egTree.tree.var #VAR output to be used as input for simultation by mutforge
      # example.gap.bed #masked out regions in BED (.bed) format
      # simple.par #PAR file to specify library type: insert distribution, coverage, read length, ...
      # example.fna #nploid sequence fasta input, use comma to separate if multiple
    
    	### simulate clonal evolution using tree-based VAR file by mutforge
      echo "### simulate clonal evolution using tree-based VAR file by mutforge"
      PYTHONPATH=$SVENGINEPATH python2 -B -m mf.mutforge -x bam -n 16 -b 100000 -f 10000 -e 1000000 -v egTree.tree.var -s example.gap.bed example.fna egTree.par example.fna
      samtools index egTree.tree.lib1.bam
    
      ## output files:
      # egTree.tree.lib1.bam #bam output, test_mf.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
      # egTree.tree.lib1.fq1/2 #fastq output, test_mf.lib2.fq1/2, ... if have multiple libs sepcified in PAR (.par) file
      # egTree.tree.out.bed #bed output for ground truth

text plot of the example tree:

|TREE|

Example of SVEngine install (Mac with home brew):
---------------------------------------------------------

To run this example, run the 'install_mac.sh' script in the 'test' subdir. Home brew is required for this provided example, although there are other ways for installation. We recommend use this script within a virtual python environment (see *virtualenv* https://pypi.python.org/pypi/virtualenv) to separate user python for SVEngine from site python. This script is easilyadapted to Ubuntu by replacing the brew install command with 'sudo apt-get install python' and also prefix all pip install commands with 'sudo'.

$SVENGINEPATH/test/install_mac.sh

  .. code:: bash

    #!/bin/bash
    
    ## install python
    ##   outside a virtual python
    brew install python2
    pip2 install virtualenv #optional
    ## install python dependencies
    ##   optionally, inside a virtual python
    pip2 install -U scipy
    pip2 install -U pysam
    pip2 install -U biopython
    pip2 install -U pybedtools
    pip2 install -U pandas
    pip2 install -U dendropy
    pip2 install -U cython
    pip2 install -U tables
    ## install svengine
    git clone https://bitbucket.org/charade/svengine.git
    git submodule update --init --recursive
    cd svengine
    python2 setup.py install

Examples of spiked-in events:
----------------------------------------------------------- 

	Showing in each figure are plottings of four sequencing feature tracks from a simulated BAM file : **read coverage**, **read pair insert size**, **hanging read pairs** and **soft clipped reads** of the region harboring the spiked-in SV. X-axis is genomic position. The background genome is simulated at 30x coverage. Two dashed virtical bar were used to delinate the SV region boundry (+/- 200bp).	 

Example Deletion (DEL):

.. |DEL| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/DEL.png
   :alt: DEL.png
   :width: 200px

|DEL|

Example Duplication (DUP):

.. |DUP| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/DUP.png
   :alt: DUP.png
   :width: 200px

|DUP|

Example Domestic Insertion (DINS):

.. |DINS| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/DINS.png
   :alt: DINS.png
   :width: 200px

|DINS|

Example Foreign Insertion (FINS):

.. |FINS| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/FINS.png
   :alt: FINS.png
   :width: 200px

|FINS|

Example Inversion (INV):

.. |INV| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/INV.png
   :alt: INV.png
   :width: 200px

|INV|

Example Inverted Domestic Insersion (IDINS):

.. |IDINS| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/IDINS.png
   :alt: IDINS.png
   :width: 200px

|IDINS|

Example Inverted Duplication (IDUP):

.. |IDUP| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/IDUP.png
   :alt: IDUP.png
   :width: 200px

|IDUP|

Example Translocation (TRA):

.. |TRAD| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/TRAD.png
   :alt: TRAD.png
   :width: 200px
.. |TRAI| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/TRAI.png
   :alt: TRAI.png
   :width: 200px

|TRAD|
|TRAI|

Example Inverted Translocation (ITRA):

.. |ITRAD| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/ITRAD.png
   :alt: ITRAD.png
   :width: 200px
.. |ITRAI| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/ITRAI.png
   :alt: ITRAI.png
   :width: 200px

|ITRAD|
|ITRAI|
