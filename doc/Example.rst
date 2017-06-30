.. contents:: Navigator

Examples
========

To easy the initial learning, we include multiple BASH scripts of SVEngine .
The scripts represent common analysis tasks with prepared input files. 

To run following examples:

  - Download and unzip SVEngine package or clone the source repository to $SVENGINEPATH.
  - Download and unzip the full test data set `svengine.test.full.tgz: <https://s3-us-west-2.amazonaws.com/lixiabucket/svengine.test.full.tgz>`__ and place the data into $SVENGINEPATH/test
  - Export the $SVENGINEPATH to your SHELL: 'export SVENGINEPATH=$SVENGINEPATH'
  - Change directory into $SVENGINEPATH/test and run the test_dep.sh to make sure dependencies OK
  - run test_fa.sh, test_mf.sh and test_tv.sh scripts for example tasks
  - required input and output files were explained inline withe scripts and in this wiki
  - run clean.sh to clear outputs from previous runs

An example SVEngine fasforge run 
--------------------------------------------------
To run this example, run the 'test_fa.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir, or it can be downloaded at: 
`svengine.test.full.tgz: <https://s3-us-west-2.amazonaws.com/lixiabucket/svengine.test.full.tgz>`__

$SVENGINEPATH/test/test_fa.sh

  .. code:: bash

    ### A simple sanity checking script for fasforge ###
    
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
      # example.fna #fasta input 
      # test_fa.meta #META file to specify event distribution
      # test_fa.var #VAR file to specify exact events
    
    	## simulating using a META file (.meta) input
      PYTHONPATH=$SVENGINEPATH python -m mf.fasforge -n 16 -b 10000 -f 10000 -e 10000 -m test_fa.meta example.gap.bed simple.par example.fna
      PYTHONPATH=$SVENGINEPATH python -m mf.fasforge --outbam -n 16 -b 10000 -f 10000 -e 10000 -m test_fa.meta example.gap.bed simple.par example.fna
    
    	## simulating using a VAR file (.var) input
      PYTHONPATH=$SVENGINEPATH python -m mf.fasforge -n 1 -b 10000 -f 10000 -e 10000 -v test_fa.var example.gap.bed simple.par example.fna
      PYTHONPATH=$SVENGINEPATH python -m mf.fasforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -v test_fa.var example.gap.bed simple.par example.fna
    
    	## simulating using both a VAR file (.var) and a META file (.meta) input (experimental)
      #PYTHONPATH=$SVENGINEPATH python -m mf.fasforge -n 16 -b 100 -f 100 -m test_fa.meta -v test_fa.var example.gap.bed simple.par example.fna
    
    	## output files:
      # test_fa.lib1.bam #bam output, test_fa.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
      # test_fa.out.bed #bed output for ground truth 
      # test_fa.out.var #var output for ground truth 
      # test_fa.fnv #fasta output after ground truth spiked-in 

An example SVEngine mutforge run 
--------------------------------------------------

To run this example, run the 'test_mf.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir, or it can be downloaded at: 
`svengine.test.full.tgz: <https://s3-us-west-2.amazonaws.com/lixiabucket/svengine.test.full.tgz>`__

$SVENGINEPATH/test/test_mf.sh

  .. code:: bash

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
    
    	## simulating using a META file (.meta) input
      PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 10000 -f 10000 -e 10000 -m test_mf.meta -s example.gap.bed example.fna simple.par example.fna
      PYTHONPATH=$SVENGINEPATH python -m mf.mutforge --outbam -n 16 -b 10000 -f 10000 -e 10000 -m test_mf.meta -s example.gap.bed example.fna simple.par example.fna
    
    	## simulating using a VAR file (.var) input
      PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 1 -b 10000 -f 10000 -e 10000 -v test_mf.var -s example.gap.bed example.fna simple.par example.fna
      PYTHONPATH=$SVENGINEPATH python -m mf.mutforge --outbam -n 1 -b 10000 -f 10000 -e 10000 -v test_mf.var -s example.gap.bed example.fna simple.par example.fna
    
    	## simulating using both a VAR file (.var) and a META file (.meta) input (experimental)
      #PYTHONPATH=$SVENGINEPATH python -m mf.mutforge -n 16 -b 100 -f 100 -m test_mf.meta -v test_mf.var example.fna example.gap.bed simple.par example.fna
    
    	## output files:
      # test_mf.lib1.bam #bam output, test_mf.lib2.bam, ... if have multiple libs sepcified in PAR (.par) file
      # test_mf.out.bed #bed output for ground truth 
      # test_mf.out.var #var output for ground truth 


An example SVEngine tree2var run 
--------------------------------------------------

To run this example, run the 'test_tv.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir, or it can be downloaded at: 
`svengine.test.full.tgz: <https://s3-us-west-2.amazonaws.com/lixiabucket/svengine.test.full.tgz>`__

$SVENGINEPATH/test/test_tv.sh

  .. code:: bash

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

text plot of the example tree:

                                                                   /------------------------------- C1
    /--------------------------------------------------------------V5                                 
    |                                                              \------------------------------- C2
    |                                                                                                 
    V1                                                             /------------------------------- C3
    |                              /-------------------------------V4                                 
    |                              |                               \------------------------------- C4
    \------------------------------V2                                                                 
                                   |                               /------------------------------- C5
                                   \-------------------------------V3                                 
                                                                   \------------------------------- C6
                                                                                                      
                                                                                                      


Example of SVEngine install (Mac with home brew):
---------------------------------------------------------

To run this example, run the 'install_mac.sh' script in the 'test' subdir. Home brew is required for this provided example, although there are other ways for installation. We recommend use this script within a virtual python environment (see *virtualenv* https://pypi.python.org/pypi/virtualenv) to separate user python for SVEngine from site python. This script is easilyadapted to Ubuntu by replacing the brew install command with 'sudo apt-get install python' and also prefix all pip install commands with 'sudo'.

$SVENGINEPATH/test/install_mac.sh

  .. code:: bash

    #!/bin/bash
    
    ## install python
    ##   outside a virtual python
    brew install pythona
    pip install virtualenv #optional
    ## install python dependencies
    ##   optionally, inside a virtual python
    pip install -U scipy
    pip install -U pysam
    pip install -U biopython
    pip install -U pybedtools
    pip install -U pandas
    pip install -U dendropy
    pip install -U cython
    pip install -U tables
    ## install svengine
    git clone https://bitbucket.org/charade/svengine.git
    git submodule update --init --recursive
    cd svengine
    python setup.py install

Examples of spiked-in events:
----------------------------------------------------------- 

	Showing in each figure are plottings of four sequencing feature tracks from a simulated BAM file : **read coverage**, **read pair insert size**, **hanging read pairs** and **soft clipped reads** of the region harboring the spiked-in SV. X-axis is genomic position. The background genome is simulated at 30x coverage. Two dashed virtical bar were used to delinate the SV region boundry (+/- 200bp).	 

Example Deletion (DEL):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/DEL.png
	 :width: 400px

Example Duplication (DUP):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/DUP.png
	 :width: 400px

Example Domestic Insertion (DINS):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/DINS.png
	 :width: 400px

Example Foreign Insertion (FINS):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/FINS.png
	 :width: 400px

Example Inversion (INV):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/INV.png
	 :width: 400px

Example Inverted Domestic Insersion (IDINS):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/IDINS.png
	 :width: 400px

Example Inverted Duplication (IDUP):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/IDUP.png
	 :width: 400px

Example Translocation (TRA):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/TRAD.png
	 :width: 400px
.. image:: https://bytebucket.org/charade/svengine/wiki/images/TRAI.png
	 :width: 400px

Example Inverted Translocation (ITRA):

.. image:: https://bytebucket.org/charade/svengine/wiki/images/ITRAD.png
	 :width: 400px
.. image:: https://bytebucket.org/charade/svengine/wiki/images/ITRAI.png
	 :width: 400px

