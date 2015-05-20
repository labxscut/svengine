#!/usr/bin/env python

"""MutForge Package.
"""

from setuptools import setup, find_packages
#from distutils.core import Extension
#from distutils.command import build
import os, sys, subprocess

doclines=__doc__.splitlines()

#class my_build(build.build): 
    # different order: build_ext *before* build_py, so that 
    # build_py can already use ctypes! 
#    sub_commands = [('build_ext', build.build.has_ext_modules), 
#        ('build_py', build.build.has_pure_modules), 
#        ('build_clib', build.build.has_c_libraries), 
#        ('build_scripts', build.build.has_scripts), ] 

#install xwgsim
os.chdir("xwgsim")
os.system("gcc -g -O2 -Wall -o xwgsim xwgsim.c -lz -lm")
os.system("mv xwgsim ~/scripts")
os.chdir("..")

setup(name="svengine",
    version="0.0.1",
    description=doclines[0],
    long_description="\n".join(doclines[2:]),
    author="Li Xia",
    author_email="li.xia@stanford.edu",
    url="https://bitbuket.com/svengine",
    license="TBD",
    platforms=["Linux"],
    packages=find_packages(exclude=['ez_setup', 'test', 'doc']),
    include_package_data=True,
    zip_safe=False,
    install_requires=["python >= 2.7","pysam >= 0.7","pybedtools >= 0.1","biopython >= 0.1","pandas >= 0.1","cython >= 0.22"],
    provides=['svengine'],
    #can insert C++ extension here when needed
    #ext_modules = [ Extension('lsa._compcore', 
    #                sources = ['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'],
    #                                          depends = ['lsa/compcore.hpp'],
    #               )],
    py_modules = [
	#'svengine.bs.asmregion','svengine.bs.mutableseq','svengine.bs.parseamos',
	#'svengine.bs.replacereads',
	'svengine.mf.mutforge'],
    #cmdclass = {'build': my_build},
    data_files = [('',['README','LICENSE.txt'])],
    entry_points = { 
      'console_scripts': [
	#'addindel = svengine.bs.addindel:run',
	#'addsv = svengine.bs.addsv:run',
	#'addsnv = svengine.bs.addsnv:run',
	'svengine = svengine.mf.svengine:run',
	'fasforge = svengine.mf.fasforge:run',
	'var2vcf = svengine.mf.var2vcf:run',
	'freqsim = svengine.mf.freqsim:run'
      ]
   },
)
