#!/bin/bash

#master doc: http://dl.dropbox.com/u/35182955/development_environment.html
#start from your home directory clean: ~
#assuming using BASH
#add following to ~/.bashrc
export PATH="$HOME/usr/bin:$PATH"
export PYTHONPATH="$HOME/usr/lib/python/python2.7/site-packages:$PATH"
. .bashrc
env | grep PATH # to check
#make necessary directories
mkdir ~/bin ~/setup ~/scripts
#return to ~/setup after every block installation
cd ~/setup
#Python (not necessary if you are using system Python)
export PVER=2.7.8
wget https://www.python.org/ftp/python/$PVER/Python-$PVER.tgz
tar -zxvf Python-$PVER.tgz
cd Python-$PVER
./configure --prefix=$HOME/usr --enable-shared
make && make install
#Python Packages
export PVER=15.1
wget https://pypi.python.org/packages/source/s/setuptools/setuptools-$PVER.tar.gz#md5=10407f6737d8dc37e5310e68c1f1f6ec
tar -zxvf setuptools-$PVER.tar.gz
cd setuptools-$PVER
python2 setup.py install
~/usr/bin/easy_install cython
#SVEngine
git clone https://bitbucket.org:charade/svengine.git
cd svengine
~/usr/bin/easy_install .
