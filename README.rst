README
--------

INTRODUCTION
============
    
    SVEngine(Structural Variants Engine)
    
    Currently the package works for Linux (tested with Ubuntu) and Mac (tested with Macports).
    It might also work for Windows with Cygwin (not tested).
    SVEngine documentation is on wiki and it is available:
    http://bitbucket.org/charade/svengine/wiki

INSTALL
========
    
    *Dependencies* (and that of):
    
    ::

      Python2(>=2.7)  http://www.python.org/

      Biopython http://www.biopython.org

      Pysam  http://pysam.readthedocs.org/en/latest/
      
        Pysam prerequisites:
        Samtools (>0.19) http://samtools.org/
        (Samtools comes wit Pysam, however it has to be able to compile on your system with dafault setting. otherwise expect the "csamtools.c" not found error. e.g. libcurse is required for build samtools but often missing in some Linux.)
    
      Pybedtools(>=3.0) http://pythonhosted.org/pybedtools/
      
        Pybedtools prerequisites:
        Pandas http://pandas.pydata.org
        Cython (>=2.2) http://cython.org/
        Numpy http://numpy.org
        Sqlite3 http://sqlite.org
   
      BWA http://samtools.org/
  
    There are multiple ways to install dependencies and SVEngine. 
    The virtualenv based installation and usage offer the benefits of non-root installation and python site separation. 
    
    **virtualenv is the MOST RECOMMENDED WAY for installation**

    1. If you have root, these can be easily installed using:
    
    ::
      
      > sudo apt-get install sqlite3 libcurses-dev python-dev python-pip
      > sudo pip install -U Pysam
      > sudo pip install -U Pandas
      > sudo pip install -U Cython
      > sudo pip install -U Pybedtools
      > sudo pip install -U Biopython
      > sudo python setup.py install
      
    2. For installations without root priviledge, you need to use virtualenv python or your local python installation,
       Please see https://virtualenv.pypa.io and http://dl.dropbox.com/u/35182955/Development_environment.html 
       for details to install virtualenv or your own python. Briefly, download virtualenv from: https://pypi.python.org/pypi/virtualenv
       then unzip and change into your virtualenv directory and run:
    
    ::
       >PYTHONPATH="local/python/prefix/lib/python2.x/site-packages" python setup.py install --prefix="local/python/prefix"
       >PYTHONPATH="local/python/prefix/lib/python2.x/site-packages" local/python/prefix/bin/virtualenv sve_vpy --no-site-packages
    
       Or if your system python has virtualenv, make sure your PYTHONPATH is set to empty and follow steps below:

    ::
     
      >virtualenv sve_vpy --no-site-packages
      
      Then you can activate this virtual python for installing SVEngine:
    
    ::
      
      >source sve_vpy/bin/activate

    3. Now either your are under virtualenv or if you are using your local python distribution, 
       these dependencies will be automatically setup by SVEngine installation:
    
    ::
      
      >python setup.py install --prefix="your/own/prefix" #if using your local python under your/own/prefix
      sve_vpy> python setup.py #if using your virtualenv python in sve_vpy
    
      Now the executables will be available from "sve_vpy/bin" or "your/own/prefix/bin". If not virtualenv, remember to add it to your PATH.
    
    4. Or if you prefer system python without root, the installation is not much different.
       please ask your system admin to install pysam, pandas, pybedtools and cython for you. 
       Change your PYTHONPATH to include "your/own/prefix" i.e.

    ::
    
      export PYTHONPATH="PYTHONPATH:your/own/prefixe/lib/python2.7/site-packages/"
    
      You might encounter other permission related problems, please also consult your admin.
      please consult the author's development document to setup a working environment:
      http://dl.dropbox.com/u/35182955/Development_environment.html

EXECUTABLES
===========
	
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
