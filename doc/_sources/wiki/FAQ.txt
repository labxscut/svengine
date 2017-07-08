.. contents:: Navigator

Usage Related
=============

1. Error "Unexpected file format.  Please use tab-delimited BED, GFF, or VCF. Perhaps you have non-integer starts or ends at line 1?"
-------------------------------------------------------------------------------------------------------------------------------------

    ANS:   Your gap.bed file, your haploid sequence file and/or your reference file are not consistent in format.
           1. Make sure all the sequence ids in .bed, .fasta files are consistent (e.g.) all with "chr" or all without "chr". To correct the fasta file, you need to rename the chromosome id and reindex the fasta file.
           2. Make sure the gap bed file doesn't contain gap information for extra chromosomes not found in fasta files. e.g. if you are spiking only into "chr22", remove any gaps defined for all other chromosomes.

Installation Related
====================

1. Why I had a 'gzopen error' from bwa ?
-------------------------------------------------------------------------------------------------------------------------------------

  ANS: Please update libz version to 1.2.8 or later.

2. What is the error: ' recognized as * from samtools means?
-------------------------------------------------------------------------------------------------------------------------------------

  ANS: There is one of your bam file empty.

3. Why I get empty .fq file even non-empty reads were generated from wgsim?
-------------------------------------------------------------------------------------------------------------------------------------
  ANS: This is a bug of wgsim, which is giving empty fq file even given non-empty but only a small number of generated reads. The way to avoid this is to increase your allele fraction or coverage.

4. Why I get error, "csamtools.c not found", from installing of pysam?
-------------------------------------------------------------------------------------------------------------------------------------
  ANS: To install pysam, you need to complie Samtools by default. Please fix samtools's dependencies first. For instance, uncurse\_ may be missed in Linux system. 

5. Why I get error, "Python.h not found", from installing of cython?
-------------------------------------------------------------------------------------------------------------------------------------
  ANS: need to install python-dev on Linux. e.g. on Ubuntu try: sudo apt-get install python-dev 
