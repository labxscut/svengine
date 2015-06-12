README
--------

INTRODUCTION
========
    
    SVEngine(Structural Variants Engine)
    
    Currently the package works for Linux (tested with Ubuntu) and Mac (tested with Macports).
    It might also work for Windows with Cygwin (not tested).
    SVEngine documentation is on wiki and it is available:
    http://bitbucket.org/charade/svengine/wiki

INSTALL
========
    
    Dependencies (and that of):
    
    ::
    
    Python2(>=2.7)  http://www.python.org/
    
    Pysam  http://pysam.readthedocs.org/en/latest/
      ::
      
      Pysam prerequisites:
      Samtools (>0.19) http://samtools.org/
      (Samtools comes wit Pysam, however it has to be able to compile on your system with dafault setting. otherwise expect the "csamtools.c" not found error. e.g. libcurse is required for build samtools but often missing in some Linux.)
    
    Pybedtools(>=3.0) http://pythonhosted.org/pybedtools/
      ::
      
      Pybedtools prerequisites:
      Pandas http://pandas.pydata.org
      Cython (>=2.2) http://cython.org/
      Numpy http://numpy.org
   
    BWA http://samtools.org/
  
    If you have root, these can be easily installed using:
    
    ::
      
      > sudo pip install -U Pandas
      > sudo pip install -U Pysam
      > sudo pip install -U Pybedtools
      > sudo pip install -U Cython
      > sudo python setup.py install
      
    Or if You are using your local python distribution, these dependencies shall be automatically setup during SVEngine installation in:
    
    ::
      
      python setup.py install --prefix="your/own/prefix"
    
    Now the executables will be available from "your/own/prefix/bin", remember to add it to your PATH.
    
    Or if you are using system python without root, please ask your system admin to install pysam, pandas, pybedtools and cython for you.
    Change your PYTHONPATH to include "your/own/prefix" i.e.
    ::
    
      export PYTHONPATH="PYTHONPATH:your/own/prefixe/lib/python2.7/site-packages/"
    
    please consult the author's development document to setup a working environment:
    http://dl.dropbox.com/u/35182955/Development_environment.html

EXECUTABLES
========
	
    xwgsim				--  enhanced wgsim
    
    mutforge            --  engineer structure variants (sophisticated)
    
    fasforge            --  engineer structure variants (straightforward)

USAGE
========
    
    (1) By default all above executables will be available from ~/scripts .
    Use '-h' to read script-wise usage.

    (2) Do a Sanity check for installation
    ::
        
        cd test
        ./test_mf.sh
        ./test_fa.sh
    
WIKI
========
    http://bitbucket.org/charade/svengine/wiki/Home
    
FAQ
========
    http://bitbucket.org/charade/svengine/wiki/FAQ
    
BUG
========
    https://bitbucket.org/charade/svengine/issues

CONTACT
========
    lixia at stanford dot edu