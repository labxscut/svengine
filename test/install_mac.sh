#!/bin/bash

## install python
##   outside a virtual python
brew install pythona
pip install virtualenv #optional
## install python dependencies
##   optionally, inside a virtual python
pip install scipy
pip install pysam
pip install biopython
pip install pybedtools
pip install pandas
pip install dendropy
pip install -U cython
## install svengine
git clone https://bitbucket.org/charade/svengine.git
git submodule update --init --recursive
cd svengine
python setup.py install
