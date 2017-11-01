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
  - Download and unzip the full test data set `svengine.test.full.tgz: <REPLACE_WITH_DOWNLOADLINK>`__ and place the data into $SVENGINEPATH/test
  - Export the $SVENGINEPATH to your SHELL: 'export SVENGINEPATH=$SVENGINEPATH'
  - Change directory into $SVENGINEPATH/test and run the test_dep.sh to make sure dependencies OK
  - run test_mf.sh and test_tv.sh scripts for example tasks
  - required input and output files were explained inline withe scripts and in this wiki
  - run clean.sh to clear outputs from previous runs

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

|TREE|

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
