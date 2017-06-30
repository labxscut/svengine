.. contents:: Navigator

Examples
========

To easy the initial learning, we include multiple BASH scripts of SVEngine .
The scripts represent common analysis tasks with prepared input files. 

To run following examples:

  - Download and unzip SVEngine package or clone the source repository to $SVENGINEPATH.
  - Download and unzip the full test data set `svengine.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__ and place the data into $SVENGINEPATH/test
  - Export the $SVENGINEPATH to your SHELL: 'export SVENGINEPATH=$SVENGINEPATH'
  - Change directory into $SVENGINEPATH/test and run the test_dep.sh to make sure dependencies OK
  - run test_fa.sh, test_mf.sh and test_tv.sh scripts for example tasks
  - required input and output files were explained inline withe scripts and in this wiki
  - run clean.sh to clear outputs from previous runs

An example SVEngine fasforge run 
--------------------------------------------------
To run this example, run the 'test_fa.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir, or it can be downloaded at: 
`svengine.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

$SVENGINEPATH/test/test_fa.sh

  .. code:: bash

REPLACE_WITH_SVENGINETESTFA

An example SVEngine mutforge run 
--------------------------------------------------

To run this example, run the 'test_mf.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir, or it can be downloaded at: 
`svengine.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

$SVENGINEPATH/test/test_mf.sh

  .. code:: bash

REPLACE_WITH_SVENGINETESTMF


An example SVEngine tree2var run 
--------------------------------------------------

To run this example, run the 'test_tv.sh' script in the 'test' subdir. Use 'clean.sh' to clean up the previous outputs first. A small test data set is required and provided within the test subdir, or it can be downloaded at: 
`svengine.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__

$SVENGINEPATH/test/test_tv.sh

  .. code:: bash

REPLACE_WITH_SVENGINETESTTV

text plot of the example tree:

REPLACE_WITH_SVENGINETREEPLOT


Example of SVEngine install (Mac with home brew):
---------------------------------------------------------

To run this example, run the 'install_mac.sh' script in the 'test' subdir. Home brew is required for this provided example, although there are other ways for installation. We recommend use this script within a virtual python environment (see *virtualenv* https://pypi.python.org/pypi/virtualenv) to separate user python for SVEngine from site python. This script is easilyadapted to Ubuntu by replacing the brew install command with 'sudo apt-get install python' and also prefix all pip install commands with 'sudo'.

$SVENGINEPATH/test/install_mac.sh

  .. code:: bash

REPLACE_WITH_SVENGINEINSTMAC

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

