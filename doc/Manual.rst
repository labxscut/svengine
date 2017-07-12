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
                            digits random): /Users/charlie/svetmp_Wcfy
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

    #METATYPE:  acronym used to specify a structural variant class:
    #           Acronyms: Deletion(DEL), Insertion(INS), Inversion(INV), Tandem Duplication(DUP), Translocation(TRA)
    #           Acronyms: Domestic INS(DINS), Foreign INS(FINS), Inverted TRA(ITRA), Inverted DUP(IDUP)
    #           the same variant class name can appear in multiple rows with different configuration
    #COUNT:     total number of the variant to be spiked-in of the varaint class defined in this row
    #SIZEDIST/SEQ: length of the sequence to be deleted/inserted if domestic insertion.
    #              Or a filename containing contigs to be randomly sampled for foreign insertion.
    #              The length could be equal prob multinormial multi(size1, size2, ...) : fix_size1_size2_...
    #              one can specify any size multiple times to increase its proportion.
    #              Or other supported distribution with parameters: dist_par1_par2_...
    #              e.g. nori_mean_sd: normal distribution rounded to integer
    #                   norf_mean_sd: normal distribution
    #                   unii_low_up: uniform distribution of integer from low to up
    #                   unif_low_up: uniform distribution
    #                   expi_rate: exponential distribution rounded to integer
    #                   expf_rate: exponential distribution
    #COPYDIST:  numpy of times that the insertion sequence to be duplicated
    #           supports applicable distribution as string encoding as SIZEDIST field
    #DONORDIST: a donor is the segment of haploid content marked and removed to be changed by the variant
    #           e.g. a deletion donates a segment without a reception
    #           e.g. a inversion donates a segment and receives a reverse complemented segment
    #           e.g. a tandem duplication donates a segment and receives a duplicted segment
    #           e.g. a insertion only receives either a domestic or a foreign segment without donation 
    #           potential values: 0, 1, 2, 3, ..., 2^nploid-1; 
    #           index to specify the group of haploids impacted by donation action; 0 is none; 2^nploid-1 is all
    #           if a sequqnce DONOR is multiple haploid then a random one among them is chosen 
    #           supports applicable distribution as string encoding as SIZEDIST field
    #RECPTDIST: a recepient is the segment of haploid receiving pontentially altered donor content
    #           potential values: 0, 1, 2, 3, ..., 2^nploid-1; 
    #           index to specify the group of haploids impacted by reception action; 0 is none; 2^nploid-1 is all
    #           if a sequence RECPT is multiple then every in the multiple is receiving the chosen donor seq
    #           e.g. FOR HOMOZYGOUS EVENTS at nploid=2: 
    #                DONOR must be fix_3, RECPT must be fix_3
    #                FOR HEMIZYGOUS EVENTS at nploid=2: 
    #                DONOR can be either fix_1 or fix_2, RECPT can be either fix_1 or fix_2
    #           supports applicable distribution as string encoding as SIZEDIST field
    #DONORFREQ: allelic frequence for donor actions (the content removing part of a structural variation) 
    #RECPTREQ:  allelic frequence for recepient actions (the content insertion part of a structural variation) 
    #NOTE:      any human readable notes go here, which will be ignored by machine
    #METATYPE	COUNT	SIZEDIST/SEQ	COPYDIST	DONORDIST	RECPTDIST	DONORFREQ	RECPTFREQ	NOTE
    FINS	2	example.fas	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#FINS=COPY*INS, foreign insertion
    ITRA	2	fix_1000	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#ITRA=1*DEL + COPY * INS, inverted transposition
    DEL	4	fix_500_500_1000_2000_2000	fix_1	fix_3	fix_3	fix_1.0	fix_0.5	#DEL=1*DEL, deletion
    TRA	2	fix_1000	fix_1	fix_3	fix_3	fix_1.0	fix_0.5	#TRA=1*DEL + COPY*INS, transposition
    INV	2	fix_1000	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#INV=1*DEL + COPY*INS, inversion
    IDUP	2	fix_1000	fix_2	fix_3	fix_3	fix_0.5	fix_0.5	#IDUP=1*DEL + COPY*INS, inverted duplication
    DINS	2	fix_1000	fix_1	fix_3	fix_3	fix_0.5	fix_0.5	#DINS=COPY*INS, domestic insertion
    DUP	2	fix_1000	fix_2	fix_3	fix_3	fix_0.5	fix_0.5	#DUP=COPY*INS, duplication

$SVENGINEPKG/svengine/test/example.par.

.. code:: bash

    #example PAR file in python config format
    #currently xwgsim is the only configurable section
    #see https://wiki.python.org/moin/ConfigParserExamples for python configure format examples
    [xwgsim]
    #number of sequencing libraries
    nlibrary=2
    #library-wise names, json format 
    libnames=["lib1","lib2"]
    #library-wise coverage, json format
    coverage: {"lib1":40,"lib2":40}
    #library-wise insert-size (isize), allowing mixtures
    #  e.g. Here we have two libraries: 
    #  lib1 has two components: [ [mode11_weight=70%, mode11_mean=250, mode11_sd=20], 
    #                             [mode12_weight=30%, mode12_mean=400, mode12_sd=40] ]
    #  lib2 has one component:  [ [mode21_weight=100%, mode21_mean=300, mode1_sd=30] ]
    isize: {"lib1":[[0.7,250,20],[0.3,400,40]],"lib2":[[1,300,30]]}
    #libarary-wise read1 length 
    read1: {"lib1":100,"lib2":100}
    #libarary-wise read2 length 
    read2: {"lib1":100,"lib2":100}
    #additional parameter to xwgsim, will be joined into a string spaced by blank
    #  use xwgsim -h to view available options
    par: {"-e":0.01,"-r":0,"-R":0,"-X":0}

$SVENGINEPKG/svengine/test/egTree.newick.

The tree input to SVEngine is in standard NEWICK format  – a widely accepted format using parenthesis to encode nested tree structures. Each internal node is labelled by the population splitting variant and weighted by the conditional splitting fraction. Each leaf node is labelled by associated terminal genotype and weighted by the subpopulation fraction (optional). See http://evolution.genetics.washington.edu/phylip/newicktree.html. For instance, the NEWICK string for the example binary tree in SVEngine's publication is: 

.. code:: bash

    ((C1,C2) V5:0.8, ((C3,C4)V4:0.8, (C5,C6)V3:0.9) V2:0.6) V1:0.5;

$SVENGINEPKG/svengine/test/egTree.newick.

.. code:: bash

    #The basic VAR file format is explained in example.var
    #The VARFREQ field in this file is dummy and to be replaced after tree fraction computation
    #It is important to note that the MID of this file must corresponds to the internal node label in the egTree.newick file
    #The MID column is used as an identifier by tree2var for correctly substituting dummy fractions
    #VID	MID	VARFREQ	VARHAP	VARCHR	VARPOS	VARDEL	VARDELSIZE	VARINS	VARINSSEQ(HAP/SEQFILE,CHR:START-END,COPY,REVCOMP)
    VAR_1	V1	1.0	0	2	1200001	True	50000	False	None
    VAR_2	V2	1.0	0	2	1300001	True	50000	False	None
    VAR_3	V3	1.0	0	2	1400001	True	50000	False	None
    VAR_4	V4	1.0	0	2	1500001	True	50000	False	None
    VAR_5	V5	1.0	0	2	1600001	True	50000	False	None
