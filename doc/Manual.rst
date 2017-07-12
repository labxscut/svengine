.. contents:: Navigator

Manual
========

0. $SVENGINEPKG environmental variable
------------------------------------------------------------------------

$SVENGINEPKG should be set to the path where you unzip or cloned the svengine package.
Below is where you can locate the original script or example file within $SVENGINEPKG.

1. Short manuals for Python and R scripts are also available with '-h' option
-------------------------------------------------------------------------

$SVENGINEPKG/svengine/mf/mutforge.py. 

The command implements a parallelized algorithm that divides template genome into trunks of contigs, spikes structural variants into contigs, samples short reads from the altered contigs, and finally merges back the short-read sets into one and performs the alignment. 

.. code:: bash

    usage: mutforge.py [-h] [-m METAFILE] [-v VARFILE] [-b BURNIN] [-e TRUNKSIZE]
                       [-o OPREFIX] [-y NPLOIDY] [-f PLANSIZE] [-g NLIGATION]
                       [-n NPROCS] [-r MERGEMAX] [-d TMPDIR] [-s SKIPFILE]
                       [-t TARGETFILE] [-x {fasta,fastq,bam}] [-c CHOOSESEED]
                       [--layout] [--debug]
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
                            var region in simulation, default: 10000000 bp
      -o OPREFIX, --oprefix OPREFIX
                            prefix for output files, default: $varfile or
                            $metafile
      -y NPLOIDY, --nploidy NPLOIDY
                            genome ploidity default= 2 ploid
      -f PLANSIZE, --plansize PLANSIZE
                            plansize + ligation*2 is the grain size for variant
                            planning, recommend plansize <
                            (genomesize/numer_of_variant/100), default: 100000 bp)
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
                            digits random): /Users/charlie/svetmp_Tjdv
      -s SKIPFILE, --skipfile SKIPFILE
                            inputfile of regions to skip, e.g. gap regions in bed
                            format
      -t TARGETFILE, --targetfile TARGETFILE
                            inputfile of regions to target, for example Exome
                            region in bed format
      -x {fasta,fastq,bam}, --chooseOutput {fasta,fastq,bam}
                            choose the type of output of mutforge, (default:
                            fastq)
      -c CHOOSESEED, --chooseSeed CHOOSESEED
                            choose the seed for python random function, (default:
                            None bp)
      --layout              dry run to layout
      --debug               debug run

$SVENGINEPKG/svengine/mf/tree2var.py. 

The tree2var command implements a procedure that determines variant fractions from an input phylogeny tree based on Equations (1) and (2) and a depth first search graph algorithm, and then substitutes these allele fractions in an input VAR file. (see SVEngine publication).

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

2. Short manuals for text based input are also found within the commented lines of the file
-------------------------------------------------------------------------------------------

$SVENGINEPKG/svengine/test/example.var.
The VAR format is intended for specifying exact information for individual variants, which includes variant id, parent id (if it is part of a complex event, e.g., the deletion part within a translocation), fraction, ploidy, chromosome, starting position, and the sequence length to be deleted and/or the sequence content to be inserted.

.. code:: bash

    #This example VAR file specified some prototype structural variants shown in the example feature track plots
    #Acronyms: Deletion(DEL), Insertion(INS), Inversion(INV), Tandem Duplication(DUP), Translocation(TRA)
    #Acronyms: Domestic INS(DINS), Foreign INS(FINS), Inverted TRA(ITRA), Inverted DUP(IDUP), 
    #VID: Variant ID for haploid-wise element variants, one unique id for each row 
    #MID: Mega ID for mega variants, e.g. a translocation with have one MID and two VIDs for del and ins parts.
    #VARFREQ: haploid-wise variant allele frequency
    #VARHAP: variant affect which haploid, 0-based, 0=first haploid 
    #VARCHR: variant affect which chromosome
    #VARDEL: if the variant carrying a deletion from its starting position, True/False only
    #VARDELSIZE: the size of the carried deletion. if no deletion, always 0 
    #VARINS: if the variant carrying a insertion from its starting position, True/False only
    #VARINSSEQ: the content of the carried insertion. if no insertion, always None.
    #           (HAP/SEQFILE,CHR:START-END,COPY,REVCOMP) - format to specify insertion content, e.g.
    #           example.fas,fins1:1-1200,1,f: says the sequence is from '1-1200bp' of 'fins1' contig in 'example.fas'.
    #                                         it is duplicated 1 times while it is not reverse complemented.
    #VID    MID     VARFREQ VARHAP  VARCHR  VARPOS  VARDEL  VARDELSIZE      VARINS  VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)
    VAR_1   FINS_1  0.5     0       1       1038000 False   0       True    example.fas,fins1:1-1200,1,f
    VAR_2   FINS_1  0.5     1       1       1038000 False   0       True    example.fas,fins1:1-1200,1,f
    VAR_3   FINS_2  0.5     0       1       1142000 False   0       True    example.fas,fins2:1-1200,1,f
    VAR_4   FINS_2  0.5     1       1       1142000 False   0       True    example.fas,fins2:1-1200,1,f
    VAR_5   ITRA_1  0.5     0       2       102000  True    1000    False   None
    VAR_6   ITRA_1  0.5     1       2       102000  True    1000    False   None
    VAR_7   ITRA_1  0.5     0       1       1246000 False   0       True    example.fna,2:102001-103000,1,r
    VAR_8   ITRA_1  0.5     1       1       1246000 False   0       True    example.fna,2:102001-103000,1,r
    VAR_9   ITRA_2  0.5     0       1       1350000 True    1000    False   None
    VAR_10  ITRA_2  0.5     1       1       1350000 True    1000    False   None
    VAR_11  ITRA_2  0.5     0       2       206000  False   0       True    example.fna,1:1350001-1351000,1,r
    VAR_12  ITRA_2  0.5     1       2       206000  False   0       True    example.fna,1:1350001-1351000,1,r
    VAR_13  DEL_1   1.0     0       1       1454000 True    1000    False   None
    VAR_14  DEL_1   1.0     1       1       1454000 True    1000    False   None
    VAR_15  DEL_2   1.0     0       2       310000  True    2000    False   None
    VAR_16  DEL_2   1.0     1       2       310000  True    2000    False   None
    VAR_17  DEL_3   1.0     0       2       414000  True    2000    False   None
    VAR_18  DEL_3   1.0     1       2       414000  True    2000    False   None
    VAR_19  DEL_4   1.0     0       2       518000  True    500     False   None
    VAR_20  DEL_4   1.0     1       2       518000  True    500     False   None
    VAR_21  TRA_1   1.0     0       2       622000  True    1000    False   None
    VAR_22  TRA_1   1.0     1       2       622000  True    1000    False   None
    VAR_23  TRA_1   0.5     0       2       726000  False   0       True    example.fna,2:622001-623000,1,f
    VAR_24  TRA_1   0.5     1       2       726000  False   0       True    example.fna,2:622001-623000,1,f
    VAR_25  TRA_2   1.0     0       2       830000  True    1000    False   None
    VAR_26  TRA_2   1.0     1       2       830000  True    1000    False   None
    VAR_27  TRA_2   0.5     0       2       934000  False   0       True    example.fna,2:830001-831000,1,f
    VAR_28  TRA_2   0.5     1       2       934000  False   0       True    example.fna,2:830001-831000,1,f
    VAR_29  INV_1   0.5     0       2       1038000 True    1000    True    example.fna,2:1038001-1039000,1,r
    VAR_30  INV_1   0.5     1       2       1038000 True    1000    True    example.fna,2:1038001-1039000,1,r
    VAR_31  INV_2   0.5     0       2       1142000 True    1000    True    example.fna,2:1142001-1143000,1,r
    VAR_32  INV_2   0.5     1       2       1142000 True    1000    True    example.fna,2:1142001-1143000,1,r
    VAR_33  IDUP_1  0.5     0       2       1246000 True    1000    True    example.fna,2:1246001-1247000,2,r
    VAR_34  IDUP_1  0.5     1       2       1246000 True    1000    True    example.fna,2:1246001-1247000,2,r
    VAR_35  IDUP_2  0.5     0       2       1350000 True    1000    True    example.fna,2:1350001-1351000,2,r
    VAR_36  IDUP_2  0.5     1       2       1350000 True    1000    True    example.fna,2:1350001-1351000,2,r
    VAR_37  DINS_1  0.5     0       2       1558000 False   0       True    example.fna,2:1454001-1455000,1,f
    VAR_38  DINS_1  0.5     1       2       1558000 False   0       True    example.fna,2:1454001-1455000,1,f
    VAR_39  DINS_2  0.5     0       2       1766000 False   0       True    example.fna,2:1662001-1663000,1,f
    VAR_40  DINS_2  0.5     1       2       1766000 False   0       True    example.fna,2:1662001-1663000,1,f
    VAR_41  DUP_1   0.5     0       2       1870000 True    1000    True    example.fna,2:1870001-1871000,2,f
    VAR_42  DUP_1   0.5     1       2       1870000 True    1000    True    example.fna,2:1870001-1871000,2,f
    VAR_43  DUP_2   0.5     0       2       1974000 True    1000    True    example.fna,2:1974001-1975000,2,f
    VAR_44  DUP_2   0.5     1       2       1974000 True    1000    True    example.fna,2:1974001-1975000,2,f


$SVENGINEPKG/svengine/test/example.meta.
The META format is intended for an easy higher-level control, allowing one to specify a desired meta distribution of variants, including variant type, and total number of events, size, allele fraction, and ploidy distributions per type. It is also possible to specify where and how to sample insertion sequence in the case of foreign insertion. For example, a user can easily request 100 deletions of size ranging from 100bp to 10kbp, of uniform allelic fraction distribution and fair coin Bernoulli distribution of homo- and heterozygosity by one line of test in the META file.  These parameters can be stipulated in a single line.

.. code:: bash

    #testing meta table
    #DONORDIST 0, 1, 2, 3, ..., 2^nploid-1; this is basically deletion or nothing
    #RECPTDIST 0, 1, 2, 3, ..., 2^nploid-1; this is basically insertion or nothing
    #IF DONOR is multiple and is providing sequence, a random one among them is chosen 
    #IF RECPT is multiple and every in the multiple is receiving the chosen donor seq
    #FOR HOMOZYGOUS EVENTS TO HAPPEN at nploid=2: 
    #DONOR must be fix_3, RECPT must be fix_3
    #And Sequence Content Should be the same for any donor sequence
    #FOR HEMIZYGOUS EVENTS TO HAPPEN at nploid=2: 
    #DONOR must be fix_1/fix_2, RECPT must be fix_1/fix_2
    #And Sequence Content Should be the same for any donor sequence
    #METATYPE	COUNT	SIZEDIST/SEQ	COPYDIST	DONORDIST	RECPTDIST	DONORFREQ	RECPTFREQ	NOTE
    FINS	2	example.fas	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#FINS=COPY*INS, foreign insertion
    ITRA	2	fix_1000	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#ITRA=1*DEL + COPY * INS, inverted transposition
    DEL	4	fix_500_500_1000_2000_2000	fix_1	fix_3	fix_3	fix_1.0	fix_0.5	#DEL=1*DEL, deletion
    TRA	2	fix_1000	fix_1	fix_3	fix_3	fix_1.0	fix_0.5	#TRA=1*DEL + COPY*INS, transposition
    INV	2	fix_1000	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#INV=1*DEL + COPY*INS, inversion
    IDUP	2	fix_1000	fix_2	fix_3	fix_3	fix_0.5	fix_0.5	#IDUP=1*DEL + COPY*INS, inverted duplication
    DINS	2	fix_1000	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#DINS=COPY*INS, domestic insertion
    DUP	2	fix_1000	fix_2	fix_3	fix_3	fix_0.5	fix_0.5	#DUP=COPY*INS, duplication

$SVENGINEPKG/svengine/test/simple.par.

