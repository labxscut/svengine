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
      

    *Installation*
    
    There are multiple ways to install dependencies and SVEngine. 
    
    The virtualenv based installation and usage offer the benefits of non-root installation and python site separation. 
    
    (1) Install SVEngine through system/site python and virtualenv
      
      **This is the MOST RECOMMENDED WAY for installation**
      
      (1.1) virtualenv command is standard with Python 2.7 or later. If it is not present, please see https://virtualenv.pypa.io for details to install virtualenv for your python. Possibly as simple as:
      
      ::

        sudo easy_install pip
        sudo pip install virtualenv

      Ask your IT manager to help install it for you if you have permission difficulties.
      
      (1.2) When your system python has virtualenv, make sure your $PYTHONPATH is set to empty and follow steps below:
      
      ::
        
        >virtualenv sve_vpy --no-site-packages

      (1.3) Then you can activate this virtual python for installing SVEngine:
      
      ::
        
        >source sve_vpy/bin/activate

      (1.4) Now under your virtualenv, the dependencies will be automatically setup by SVEngine installation itself:
      
      ::
        
        sve_vpy> python setup.py install 

      (1.5) On some platforms, you might run into error and it prompts you to install Cython, just install Cython and try again:
      
      ::
        
        sve_vpy>pip install -U cython
      
      (1.6) Now the SVEngine executables will be available from "$PWD/sve_vpy/bin". Because you installed SVEngine via virtualenv, remember to activate the virtualenv first every time you use SVEngine. Also export the environmental variable $SVENGINE_BIN=$PWD/sve_vpy/bin
       
    (2) Other ways to install SVEngine

      (2.1) If you have root, these can be easily installed using:
      
      ::
        
        > sudo apt-get install sqlite3 libcurses-dev python-dev python-pip
        > sudo pip install -U Pysam
        > sudo pip install -U Pandas
        > sudo pip install -U Cython
        > sudo pip install -U Pybedtools
        > sudo pip install -U Biopython
        > sudo python setup.py install
      
      The executables are at $SVENGINE_BIN=$(dirname `which python`)
      
      (2.2) For installations without root priviledge you need your local python installation and virtualenv. Please see the Author's http://dl.dropbox.com/u/35182955/Development_environment.html for details to install your own python. Also see the install script examples in https://bitbucket.org/charade/svengine/wiki/Example for detailed steps.

      With your own python, download virtualenv from: https://pypi.python.org/pypi/virtualenv then unzip and change into your virtualenv directory and run:
      
      ::
        
        >PYTHONPATH="/local/python/prefix/lib/python2.x/site-packages" python setup.py install --prefix="/local/python/prefix"
        >PYTHONPATH="/local/python/prefix/lib/python2.x/site-packages" /local/python/prefix/bin/virtualenv sve_vpy --no-site-packages

      Then the rest steps are just exactly the same to (1.3) - (1.6).
      Or, you can elect to install SVEngine to your local python directly: 
      
      ::

        >python setup.py install --prefix="/local/python/prefix" #if using your local python under your own /local/python/prefix

      The executables are in "/local/python/prefix/bin". So set $SVENGINE_BIN="/local/python/prefix/bin" and remember to add it to your $PATH.

    3. There is also the possibility to install with your system site python without root. 
      The installation process is not much different but please first ask your system admin to install pysam, pandas, pybedtools and cython to the system site python for you. Then change your PYTHONPATH to include "/local/python/prefix" which is writable to you as your user python library path. i.e.
      
      ::

        export PYTHONPATH="PYTHONPATH:/local/python/prefix/lib/python2.7/site-packages/"
      
      The rest are the same with local python installation (either using virtualenv or not). But you might encounter other permission related problems, please consult your admin. Also, please consult the author's development document on how to setup a working environment: http://dl.dropbox.com/u/35182955/Development_environment.html

EXECUTABLES
===========

    xwgsim        --  enhanced wgsim
    
    mutforge      --  engineer structure variants 
    
    fasforge      --  engineer structure variants 


USAGE
========
    
    (1) By default all above executables will be available from ~/scripts .
    Use '-h' to read script-wise usage.

    (2) Do a Sanity check for installation

    ::

      cd test
      ./test_mf.sh
    
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
