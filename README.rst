.. |Logo| image:: https://bitbucket.org/charade/svengine/raw/master/doc/images/logo.png
   :alt: logo.png
   :height: 50px
   :width: 50px

.. |Rationale| image:: https://bitbucket.org/charade/zoomx/raw/master/doc/images/rationale.png
   :alt: rationale.png
   :height: 450px
   :width: 540px

.. |Pipeline| image:: https://bitbucket.org/charade/zoomx/raw/master/doc/images/pipeline.png
   :alt: pipeline.png
   :height: 320px
   :width: 520px

|Logo| SVEngine: Allele Specific and Haplotype Aware Structural Variants Simulator
==================================================================================

.. contents:: Navigator

QUICK LINKS
-----------

`Examples <https://bitbucket.org/charade/svengine/wiki/Example>`__

`Manuals <https://bitbucket.org/charade/svengine/wiki/Manual>`__

`FAQ <https://bitbucket.org/charade/svengine/wiki/FAQ>`__


INTRODUCTION
------------
    
    SVEngine (Structural Variants Engine)

  - SVEngine is a multi-purpose and self-contained simulator for whole genome scale spike-in of thousands of SV events of various types in both single-sample and matched sample scenarios.

  - SVEngine takes as input reference contigs in FASTA files, variant meta distribution as specified in META files (see Manual) or specific variant information as specified in VAR files (see Manual) and NEWICK files for specifying clonal phylogenetic trees in cancer.

  - SVEngine outpus alterred contigs in FASTA files, spiked-in variants in VAR files (see Manual), simulated short read in FASTQ files and aligned short reads in BAM files.

  - SVEngine's modeling and pipeline are illustrated in Figures 1 and 2.

  Figure 1. The principle and workflow of SVEngine.

  |Rationale|

  Figure 2. The bioinformatics pipeline of SVEngine.

  |Pipeline|
    
  Currently the package works for Linux (tested with Ubuntu) and Mac (tested with Homebrew).
  It might also work for Windows with Cygwin (not tested).
  SVEngine documentation is on wiki and it is available:
  http://bitbucket.org/charade/svengine/wiki

TRYOUT
------------

  An AWS EC2 instance and a preconfigured virtual machine are available for trying out SVEngine without installation efforts.

  *AWS*:

    SVEngine is available for tryout in a live Amazon EC2 instance.
    Simply log on with SSH and find the package located in ~/setup/svengine .

    .. code::

      ssh bkd@ec2-54-190-36-93.us-west-2.compute.amazonaws.com # secret is the same as user

  Testing data and scripts are available in ~/setup/svengine/test. Run ./clean.sh to clean up the directory first. Examples are offered and explained in `Examples <https://bitbucket.org/charade/svengine/wiki/Example>`__ .

  *VirtualMachine*:

    The Amazon Machine Image (ami) of the above EC2 instance is downloadable from Amazon S3 bucket to run locally.
    https://s3-us-west-2.amazonaws.com/lixiabucket/svengineDemo.ova . The usage is the same as above.

INSTALL
------------

  SVEngine is available on bitbucket: bitbucket.org:charade/svengine.git

  *Dependencies*:

    Python2(>=2.7)  http://www.python.org/

      pybedtools, pysam, pandas, tables, scipy, dendropy, cython


    A script to confirm all dependencies are met is here:

      https://bitbucket.org/charade/svengine/raw/master/test/test_dep.sh

    An example dependency install script for Mac wth Homebrew is here:

      https://bitbucket.org/charade/svengine/raw/master/test/install_mac.sh


  *Installation Free !*:

    SVEngine can be used without installation as long as required commands are all in $PATH.
    see `Examples <https://bitbucket.org/charade/svengine/wiki/Example>`__ and `Manuals <https://bitbucket.org/charade/svengine/wiki/Manual>`__ for more informations

  *Installation*

    You can install SVEngine as a standard package to your site python without virtualenv, in particular when all python dependencies are already met. Before you take this approach, please consider following virtualenv based installation first. The following guides were intended for Ubuntu, while Mac Homebrew user may easily adapt by removing sudo.

    The other option is to use virtualenv based non-root installation, which separates svengine from site python. virtualenv command is standard with Python 2.7 or later. see: http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/. If it is not present, please see https://virtualenv.pypa.io for details to install virtualenv for your python. 

    .. code::

      sudo easy_install pip
      sudo pip install virtualenv

    Ask your IT manager to help install it for you if you have permission difficulties.
    When your system python has virtualenv, make sure your $PYTHONPATH is set to empty.
    Follow steps below:

    .. code::

      sudo pip install --upgrade pip
      virtualenv svengine_vpy

    Now you can activate this virtual python for installing SVEngine python dependencies:

    .. code::

      # in svengine_vpy>
      source svengine_vpy/bin/activate
      pip install scipy
      pip install tables
      pip install pandas
      pip install pysam
      pip install pybedtools
      pip install dendropy
      pip install biopython
      pip install -U cython

    Then install svengine as an standard python package:

    .. code::

      # in svengine_vpy>
      git clone git@bitbucket.org:charade/svengine.git
      cd svengine
      git submodule update --init --recursive
      python setup.py install

    Now the SVEngine executables will be available from "$PWD/svengine_vpy/bin".

    Because you installed SVEngine via virtualenv, remember to activate the virtualenv every time you use SVEngine. Also export the environmental variable $SVENGINE_BIN=$PWD/svengine_vpy/bin and add it to your $PATH.

EXECUTABLES
-----------

    xwgsim        --  enhanced Heng Li's wgsim
    
    mutforge      --  spiking in structure variants without fasta output 
    
    fasforge      --  spiking in structure variants with fasta output  

    tree2var      --  updating variant with a phylogenetic tree structured frequency

USAGE
------------

  By default all above executables will be available from $SVENGINE_BIN/ .
  Use '-h' to read script-wise usage.
  See Wiki pages for more details.

WIKI
------------

  SVEngine's https://bitbucket.org/charade/svengine/wiki page is a growing resource for manuals, FAQs and other information. This is a MUST read place before you actually using the SVEngine tool. These documentations are also openly editable in reSructuredText format. You are more than welcome to contribute to this ongoing documentation.
    
FAQ
------------

  http://bitbucket.org/charade/svengine/wiki/FAQ
    
BUG
------------

  https://bitbucket.org/charade/svengine/issues

CONTACT
------------

  Questions and comments shall be addressed to lixia at stanford dot edu.


CITATIONS
---------

  The SVEngine manuscript is under review.
