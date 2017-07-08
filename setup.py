#!/usr/bin/env python

"""MutForge Package.
"""

from setuptools import setup, find_packages
import os, sys, subprocess, ConfigParser

doclines=__doc__.splitlines()

#install bwa, samtools, bedtools, dwgsim
#all these are bundled with SVEngine, no separate download

"""samtools"""

print "cd samtools && make"
os.system("cd samtools && make")
print "cp samtools/samtools %s" % sys.prefix
os.system("cp samtools/samtools %s/bin" % sys.prefix)

"""bwa"""

print "cd bwa && make"
os.system("cd bwa && make")
print "cp bwa/bwa %s" % sys.prefix
os.system("cp bwa/bwa %s/bin" % sys.prefix)

"""dwgsim"""

if not os.path.exists("dwgsim/samtools"):  
	print "cp -r samtools dwgsim/"
	os.system("cp -r samtools dwgsim/")
print "cd dwgsim && make"
os.system("cd dwgsim && make")
print "cp dwgsim/dwgsim %s" % sys.prefix
os.system("cp dwgsim/dwgsim %s/bin" % sys.prefix)

"""bedtools"""

print "cd bedtools && make"
os.system("cd bedtools && make")
print "cp bedtools/bin/* %s/bin" % sys.prefix
os.system("cp bedtools/bin/* %s/bin" % sys.prefix)

dist=setup(name="svengine",
    version="1.0.0",
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
    install_requires=["python >= 2.7","cython >= 0.23","pandas >= 0.16","pybedtools >= 0.6.9","pysam >= 0.7","biopython >= 0.1",],
    provides=['svengine'],
    py_modules = ['svengine.mf.mutforge'],
		scripts = ['xwgsim/xwgsim'],
    data_files = [('',['README.rst','LICENSE'])],
    entry_points = { 
      'console_scripts': [
	       'mutforge = svengine.mf.mutforge:run',
	       'fasforge = svengine.mf.fasforge:run',
	       'var2vcf = svengine.mf.var2vcf:run',
	       'tree2var = svengine.mf.tree2var:run',
	       'freqsim = svengine.mf.freqsim:run'
      ]
   },
)

#print(dist.dump_option_dicts())
#print(dist.parse_config_files())
#print(dist.metadata)

