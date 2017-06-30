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
