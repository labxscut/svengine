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
