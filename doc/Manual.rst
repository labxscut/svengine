.. contents:: Navigator

Manual
========

0. $SVENGINEPKG environmental variable
------------------------------------------------------------------------

$SVENGINEPKG should be set to the path where you unzip or cloned the svengine package.
Below is where you can locate the original script for each executables within $SVENGINEPKG

1. Short Manuals for Python and R scripts are available with '-h' option
-------------------------------------------------------------------------

$SVENGINEPKG/svengine/mf/mutforge.py

.. code:: bash

    usage: mutforge.py [-h] [-m METAFILE] [-v VARFILE] [-b BURNIN] [-e TRUNKSIZE]
                       [-i INBAMFILENAMES] [-o OPREFIX] [-y NPLOID] [-f PLANSIZE]
                       [-g NLIGATION] [-n NPROCS] [-r MERGEMAX] [-d TMPDIR]
                       [-s SKIPFILE] [-t TARGETFILE] [--layout] [--outbam]
                       [--debug]
                       hapfiles parfile reffile
    
    spike-in mutations to fasta/bam files, for bam files first assemble haplo-
    types for fasta files we segment fasta files to confirm if the configuration
    is viable
    
    positional arguments:
      hapfiles              inputfile of hap file(s), required
      parfile               inputfile of library-wise parameters, required
      reffile               inputfile of ref file(s), required
    
    optional arguments:
      -h, --help            show this help message and exit
      -m METAFILE, --metafile METAFILE
                            inputfile of meta specification
      -v VARFILE, --varfile VARFILE
                            inputfile of variant specification
      -b BURNIN, --burnin BURNIN
                            burnin is the chromosome tip skip size for genome
                            planning to avoid tip padding Ns, default: 1100000 bp
      -e TRUNKSIZE, --trunksize TRUNKSIZE
                            trunksize is the parallelized compartment size of non-
                            var region in simulation 10000000 bp
      -i INBAMFILENAMES, --bamfile INBAMFILENAMES
                            inputfile of library-wise bam file(s)
      -o OPREFIX, --oprefix OPREFIX
                            prefix for output files
      -y NPLOID, --nploid NPLOID
                            genome ploidity default= 2 ploid
      -f PLANSIZE, --plansize PLANSIZE
                            plansize + ligation*2 is the grain size for genome
                            planning, default: 100000 bp)
      -g NLIGATION, --nligation NLIGATION
                            genome trunk ligation, recommended
                            lib_is_max+2*lib_is_sd_max (default: 500 bp)
      -n NPROCS, --nprocs NPROCS
                            split into multiple processes (default: 1 )
      -r MERGEMAX, --mergemax MERGEMAX
                            max number of bamfiles in one iteration of samtool
                            merging (default: 4000 )
      -d TMPDIR, --tmpdir TMPDIR
                            root dir for keeping temporary files, default (last 4
                            digits random): /home/lixia/svetmp_yjRU
      -s SKIPFILE, --skipfile SKIPFILE
                            inputfile of regions to skip, e.g. gap regions in bed
                            format
      -t TARGETFILE, --targetfile TARGETFILE
                            inputfile of regions to target, for example Exome
                            region in bed format
      --layout              dry run to layout
      --outbam              wet run to bam
      --debug               debug run

$SVENGINEPKG/svengine/mf/tree2var.py

.. code:: bash

    usage: tree2var.py [-h] [-o OPREFIX] [-p] varfile nwkfile
    
    assign tree structures to VAR files
    
    positional arguments:
      varfile               inputfile of SVEngine's VAR format, required
      nwkfile               inputfile of Newick format where the internal nodes
                            are named and matches to MID used in varfile, required
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OPREFIX, --oprefix OPREFIX
                            prefix for output files, default is the same as
                            varfile's prefix yet change suffix to .vtr
      -p, --plottree        generate a text plot of a evolution phylogenetic tree

2. Short Manuals for Bash scripts are available with comment lines or call wihtout options  
-------------------------------------------------------------------------------------------

