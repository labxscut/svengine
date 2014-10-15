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

setup(name="mutforge",
    version="0.0.1",
    description=doclines[0],
    long_description="\n".join(doclines[2:]),
    author="Li Xia",
    author_email="li.xia@stanford.edu",
    url="https://bitbuket.com/mutforge",
    license="TBD",
    platforms=["Linux"],
    packages=find_packages(exclude=['ez_setup', 'test', 'doc']),
    include_package_data=True,
    #packages=['lsa'],
    #requires=['numpy(>=1.1)','scipy(>=0.6)','python(>=2.5)','matplotlib(>=0.98)'],
    zip_safe=False,
    install_requires=["python >= 2.7","pysam >= 0.7","pybedtools >= 0.1","biopython >= 0.1","pandas >= 0.1"],
    provides=['mutforge'],
    #can insert C++ extension when needed
    #ext_modules = [ Extension('lsa._compcore', 
    #                sources = ['lsa/compcore_wrap.cpp', 'lsa/compcore.cpp'],
    #                                          depends = ['lsa/compcore.hpp'],
    #               )],
    py_modules = ['mutforge.bs.asmregion','mutforge.bs.mutableseq','mutforge.bs.parseamos','mutforge.bs.replacereads','mutforge.mf.mutforge'],
    #cmdclass = {'build': my_build},
    data_files = [('',['README','LICENSE.txt'])],
    entry_points = { 
      'console_scripts': [
            'addindel = mutforge.bs.addindel:run',
            'addsv = mutforge.bs.addsv:run',
            'addsnv = mutforge.bs.addsnv:run',
						'mutforge = mutforge.mf.mutforge:run',
						'fasforge = mutforge.mf.fasforge:run'
      ]
   },
)
